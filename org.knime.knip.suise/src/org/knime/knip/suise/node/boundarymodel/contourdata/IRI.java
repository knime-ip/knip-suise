/*
 * ------------------------------------------------------------------------
 *
 *  Copyright by
 *  University of Konstanz, Germany and
 *  KNIME GmbH, Konstanz, Germany
 *  Website: http://www.knime.org; Email: contact@knime.org
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License, Version 3, as
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, see <http://www.gnu.org/licenses>.
 *
 *  Additional permission under GNU GPL version 3 section 7:
 *
 *  KNIME interoperates with ECLIPSE solely via ECLIPSE's plug-in APIs.
 *  Hence, KNIME and ECLIPSE are both independent programs and are not
 *  derived from each other. Should, however, the interpretation of the
 *  GNU GPL Version 3 ("License") under any applicable laws result in
 *  KNIME and ECLIPSE being a combined program, KNIME GMBH herewith grants
 *  you the additional permission to use and propagate KNIME together with
 *  ECLIPSE with only the license terms in place for ECLIPSE applying to
 *  ECLIPSE and the GNU GPL Version 3 applying for KNIME, provided the
 *  license terms of ECLIPSE themselves allow for the respective use and
 *  propagation of ECLIPSE together with KNIME.
 *
 *  Additional permission relating to nodes for KNIME that extend the Node
 *  Extension (and in particular that are based on subclasses of NodeModel,
 *  NodeDialog, and NodeView) and that only interoperate with KNIME through
 *  standard APIs ("Nodes"):
 *  Nodes are deemed to be separate and independent programs and to not be
 *  covered works.  Notwithstanding anything to the contrary in the
 *  License, the License does not apply to Nodes, you are not required to
 *  license Nodes under the License, and you are granted a license to
 *  prepare and propagate Nodes, in each case even if such Nodes are
 *  propagated with or for interoperation with KNIME.  The owner of a Node
 *  may freely choose the license terms applicable to such Node, including
 *  when such Node is propagated with or for interoperation with KNIME.
 * --------------------------------------------------------------------- *
 *
 */
package org.knime.knip.suise.node.boundarymodel.contourdata;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import net.imglib2.Pair;
import net.imglib2.util.ValuePair;

import org.knime.knip.core.util.PermutationSort;

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.AbstractClassifier;
import weka.core.Capabilities;
import weka.core.Capabilities.Capability;
import weka.core.DistanceFunction;
import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.Instances;
import weka.filters.unsupervised.attribute.MultiInstanceToPropositional;

/**
 * Interval Rule Induction
 * 
 * TODO's: contour models (consisting of different rules)
 * 
 * @author hornm, University of Konstanz
 */
public class IRI extends AbstractClassifier {

    private double m_sampleRate = .1;

    private double m_bias = 100;

    private double m_coverRate = .8;

    private List<IntervalRule> m_rules = new ArrayList<IntervalRule>();

    private int m_numThreads = 8;

    /**
     * {@inheritDoc}
     */
    @Override
    public void buildClassifier(Instances miData) throws Exception {

        // can classifier handle the data?
        getCapabilities().testWithFail(miData);

        final Instances tmpMiData = new Instances(miData);
        final Instances flatData = toSingleInstanceDataset(miData, null);

        int numPosBags = 0;
        for (Instance bag : miData) {
            if (bag.value(2) == 1) {
                numPosBags++;
            }
        }

        int remainingNumPosBags = numPosBags;
        Future<Pair<IntervalRule, Double>>[] futures = new Future[m_numThreads];
        ExecutorService pool = Executors.newFixedThreadPool(m_numThreads);

        while (remainingNumPosBags / (double)numPosBags > 1 - m_coverRate) {

            final int numIterations =
                    ((int)(m_sampleRate * remainingNumPosBags)) / m_numThreads
                            + 1;

            for (int t = 0; t < m_numThreads; t++) {
                futures[t] =
                        pool.submit(new Callable<Pair<IntervalRule, Double>>() {
                            @Override
                            public Pair<IntervalRule, Double> call()
                                    throws Exception {
                                return createRule(flatData, tmpMiData,
                                        numIterations);
                            }
                        });
            }

            // select the best rule from the threads
            double score = -Double.MAX_VALUE;
            IntervalRule rule = null;
            for (int f = 0; f < futures.length; f++) {
                if (futures[f].get().getB() > score) {
                    score = futures[f].get().getB();
                    rule = futures[f].get().getA();
                }
            }

            m_rules.add(rule);

            // only keep the bags whose instances are not covered by this rule
            Instances tmp = new Instances(tmpMiData);
            tmpMiData.clear();
            boolean covered;
            remainingNumPosBags = 0;
            for (Instance bag : tmp) {
                covered = false;
                for (Instance inst : bag.relationalValue(1)) {
                    double[] distr;
                    distr = rule.distributionForInstance(inst);
                    if (distr[1] > distr[0]) {
                        covered = true;
                        break;
                    }
                }
                if (!covered) {
                    tmpMiData.add(bag);
                    if (bag.value(2) == 1) {
                        remainingNumPosBags++;
                    }
                }
            }
            flatData.clear();
            toSingleInstanceDataset(tmpMiData, flatData);
        }

        pool.shutdown();

    }

    private Pair<IntervalRule, Double> createRule(Instances flatData,
            Instances miData, int iterations) throws Exception {
        // store for the distances between the reference distance and all others
        double[] distances = new double[flatData.numInstances()];
        // the distance function
        DistanceFunction distFunc = new EuclideanDistance(flatData);
        // permutation which sorts the distances
        Integer[] perm = new Integer[flatData.numInstances()];

        IntervalRule bestRule = null;
        double bestRuleScore = -Double.MAX_VALUE;

        // retrieve the best rule heuristically for a number of iterations
        for (int ruleIterations = 0; ruleIterations < iterations; ruleIterations++) {

            // System.out.println("------- Iteration " + ruleIterations
            // + "----------");

            // randomly select an initial instance, i.e. selecting a positive
            // bag
            // randomly and taking the instance with the largest weight
            Random r = new Random();
            int bagIdx;
            while (miData.get(bagIdx = r.nextInt(miData.numInstances())).value(
                    2) == 0)
                ;

            // the reference instance for the next rule
            Instance refInstance =
                    miData.get(bagIdx).relationalValue(1).firstInstance();
            for (Instance i : miData.get(bagIdx).relationalValue(1)) {
                if (i.weight() > refInstance.weight()) {
                    refInstance = i;
                }
            }

            // System.out.println("\tRef Instance: " + refInstance);

            IntervalRule rule = new IntervalRule();
            rule.updateClassifier(refInstance);

            // calculate the distance from that particular reference instance to
            // all other
            // positive instances (negatives are set to NaN) and sort them
            Arrays.fill(distances, Double.NaN);
            for (int i = 0; i < distances.length; i++) {
                if (flatData.get(i).classValue() == 1) {
                    distances[i] =
                            distFunc.distance(refInstance, flatData.get(i));
                }
            }
            PermutationSort.sortPermInPlace(distances, perm);

            double ruleScore = 0;
            double tmpRuleScore = 0;

            // extend the rule successively by the nearest instances till the
            // score doesn't increase anymore

            int instanceIdx = 0;
            while (true) {
                if (!Double.isNaN(distances[perm[instanceIdx]])) {
                    IntervalRule tmpRule = new IntervalRule(rule);
                    tmpRule.updateClassifier(flatData.get(perm[instanceIdx]));

                    // System.out.println("\tNext Instance: "
                    // + flatData.get(perm[instanceIdx]));
                    // System.out.println("\tCurrent Rule: " + tmpRule);

                    // evaluate rule
                    tmpRuleScore = ruleScore(tmpRule, flatData);

                    if (tmpRuleScore >= ruleScore) {
                        ruleScore = tmpRuleScore;
                        rule = tmpRule;
                    } else {
                        break;
                    }
                }
                instanceIdx++;
            }

            if (ruleScore > bestRuleScore) {
                bestRuleScore = ruleScore;
                bestRule = rule;
            }

        } // iterations per rule

        return new ValuePair<IntervalRule, Double>(bestRule, bestRuleScore);
    }

    private double ruleScore(IntervalRule rule, Instances data)
            throws Exception {

        double posCount = 0;
        double negCount = 0;
        double posSumWeights = 0;
        for (Instance inst : data) {
            if (inst.weight() > 0) {
                double dist[] = rule.distributionForInstance(inst);
                if (dist[1] > dist[0]) {
                    if (inst.classValue() == 1) {
                        posSumWeights += inst.weight();
                        posCount++;
                    } else {
                        negCount++;
                    }
                }
            }
        }
        double score = posSumWeights / (posCount + negCount + m_bias);

        // System.out.println("\tpSW=" + posSumWeights + ";pC=" + posCount
        // + ";nC=" + negCount + ";score=" + score);

        return score;
    }

    private Instances toSingleInstanceDataset(Instances miData,
            Instances flatData) throws Exception {
        MultiInstanceToPropositional convertToProp =
                new MultiInstanceToPropositional();

        convertToProp.setInputFormat(miData);

        for (int i = 0; i < miData.numInstances(); i++) {
            convertToProp.input(miData.instance(i));
        }
        convertToProp.batchFinished();

        if (flatData == null) {
            flatData = convertToProp.getOutputFormat();
            flatData.deleteAttributeAt(0); // remove the bag index attribute

        }

        Instance processed;
        while ((processed = convertToProp.output()) != null) {
            processed.setDataset(null);
            processed.deleteAttributeAt(0); // remove the bag index attribute
            flatData.add(processed);
        }

        // remove class attribute
        // flatData.setClassIndex(-1);
        // flatData.deleteAttributeAt(flatData.numAttributes() - 1);

        // set weights
        int instanceIdx = 0;
        for (Instance bag : miData) {
            for (Instance instance : bag.relationalValue(1)) {
                flatData.get(instanceIdx).setWeight(instance.weight());
                instanceIdx++;
            }
        }
        return flatData;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] distributionForInstance(Instance bag) throws Exception {
        Instances contents = bag.relationalValue(1);
        double[] res = null;
        boolean positive = false;
        for (Instance i : contents) {
            res = distributionForSingleInstance(i);
            if (res[1] > res[0]) {
                positive = true;
                break;
            }
        }
        if (positive) {
            return res;
        } else {
            return new double[]{1.0, 0};
        }
    }

    public double[] distributionForSingleInstance(Instance instance)
            throws Exception {
        double res = 0;
        for (IntervalRule r : m_rules) {
            res = Math.max(res, r.distributionForInstance(instance)[1]);
        }
        return new double[]{1 - res, res};
    }

    /**
     * @param instance
     * @return -1 if no fitting rule was found
     * @throws Exception
     */
    public int getBestRuleIndexForSingleInstance(Instance instance)
            throws Exception {
        int bestIdx = -1;
        double bestProb = 0;
        for (int i = 0; i < m_rules.size(); i++) {
            double[] distr = m_rules.get(i).distributionForInstance(instance);
            if (distr[1] > bestProb) {
                bestProb = distr[1];
                bestIdx = i;
            }
        }
        return bestIdx;

    }

    public IntervalRule getRule(int index) {
        return m_rules.get(index);
    }

    public void setBias(double bias) {
        m_bias = bias;
    }

    public void addRules(IntervalRule... rules) {
        m_rules.addAll(Arrays.asList(rules));

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Capabilities getCapabilities() {
        Capabilities result = super.getCapabilities();

        // attributes
        result.enable(Capability.NOMINAL_ATTRIBUTES);
        result.enable(Capability.RELATIONAL_ATTRIBUTES);
        result.disable(Capability.MISSING_VALUES);

        // class
        result.disableAllClasses();
        result.disableAllClassDependencies();
        result.enable(Capability.BINARY_CLASS);

        // Only multi instance data
        result.enable(Capability.ONLY_MULTIINSTANCE);

        return result;
    }

}
