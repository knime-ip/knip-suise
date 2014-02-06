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

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.Classifier;
import weka.clusterers.Clusterer;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instances;

/**
 * 
 * @author hornm, University of Konstanz
 */
public class WekaContourDataClassifier implements ContourDataClassifier {

    private Classifier m_classifier;

    //    private Instances m_data;

    private int m_numContourClusters;

    private ContourDataExtractor m_cDataSelection;

    private final Clusterer m_clusterer;

    private DenseInstance m_tmpInstance;

    private double[] m_tmpVec;

    private Instances m_data;

    /**
     * @param classifier
     */
    public WekaContourDataClassifier(Classifier classifier, ContourDataExtractor cData, Clusterer clusterer) {
        m_classifier = classifier;
        m_cDataSelection = cData;
        m_clusterer = clusterer;

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void buildClassifier(ContourDataGrid cDataGrid, VectorDataList bgData) throws Exception {
        m_cDataSelection.extractContourData(cDataGrid);
        m_tmpVec = new double[cDataGrid.numFeatures()];
        m_tmpInstance = new DenseInstance(1.0, m_tmpVec);

        if (m_clusterer != null) {
            //cluster contour data to create different models
            ArrayList<Attribute> attInfo = new ArrayList<Attribute>(m_cDataSelection.numFeatures());
            for (int i = 0; i < m_cDataSelection.numFeatures(); i++) {
                attInfo.add(new Attribute("" + i));
            }
            Instances instances = new Instances("contour_data", attInfo, m_cDataSelection.numVectors());
            for (int i = 0; i < m_cDataSelection.numVectors(); i++) {
                instances.add(new DenseInstance(m_cDataSelection.weight(i), m_cDataSelection.getVector(i)));
            }
            m_clusterer.buildClusterer(instances);

            for (int i = 0; i < m_cDataSelection.numVectors(); i++) {
                m_cDataSelection.setClusterIdx(i, m_clusterer.clusterInstance(new DenseInstance(m_cDataSelection
                        .weight(i), m_cDataSelection.getVector(i))));

            }
        }

        Instances data =
                initDataset(m_cDataSelection.numFeatures(), m_cDataSelection.numClusters() + bgData.numClusters(),
                            m_cDataSelection.numVectors() + bgData.numVectors());

        m_numContourClusters = m_cDataSelection.numClusters();

        // positive training samples
        for (int n = 0; n < m_cDataSelection.numVectors(); n++) {
            if (m_cDataSelection.weight(n) > 0) {
                addInstance(m_cDataSelection.getVector(n), m_cDataSelection.weight(n),
                            m_cDataSelection.getClusterIdx(n), data);
            } else {
                //if weight == 0, add the according instance to the negative training samples
                addInstance(m_cDataSelection.getVector(n), 1, m_cDataSelection.numClusters() + bgData.getClusterIdx(n),
                            data);
            }
        }

        // negative training samples from background
        for (int n = 0; n < bgData.numVectors(); n++) {
            if (bgData.weight(n) > 0) {
                addInstance(bgData.getVector(n), 1, m_cDataSelection.numClusters() + bgData.getClusterIdx(n), data);
            }
        }

        // negative training samples from the data grid
        for (double[] vec : m_cDataSelection.nonContourVectors()) {
            addInstance(vec, 1, m_cDataSelection.numClusters(), data);
        }

        m_classifier.buildClassifier(data);
        m_tmpInstance.setDataset(data);
        m_data = data;

    }

    private void addInstance(double[] vec, double weight, double classValue, Instances dataset) {
        DenseInstance inst = new DenseInstance(weight, vec);
        inst.insertAttributeAt(inst.numAttributes());
        inst.setDataset(dataset);
        inst.setClassValue(classValue);
        dataset.add(inst);
    }

    public ContourDataExtractor getClusteredContourData() {
        return m_cDataSelection;
    }

    /**
     * {@inheritDoc}
     * 
     * @throws Exception
     */
    @Override
    public double contourProbability(double[] inst) throws Exception {
        double[] distr = contourProbDistribution(inst);
        double max = 0;
        for (int j = 0; j < distr.length; j++) {
            max = Math.max(max, distr[j]);
        }
        return max;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] contourProbDistribution(double[] inst) throws Exception {
        //        System.arraycopy(inst, 0, m_tmpVec, 0, inst.length);
        DenseInstance dinst = new DenseInstance(1.0, inst);
        dinst.setDataset(m_data);
        double[] tmp = m_classifier.distributionForInstance(dinst);
        double[] distr = new double[m_numContourClusters];
        System.arraycopy(tmp, 0, distr, 0, distr.length);
        return distr;
    }

    /*
     * Helper method
     */
    private Instances initDataset(int numFeatures, int numClasses, int capacity) {
        ArrayList<Attribute> attributes = new ArrayList<Attribute>(numFeatures + 5);
        for (int i = 0; i < numFeatures; i++) {
            attributes.add(new Attribute("att" + i));
        }

        ArrayList<String> classNames = new ArrayList<String>();
        for (int i = 0; i < numClasses; i++) {
            classNames.add("class" + i);
        }

        Attribute classAtt = new Attribute("class", classNames);
        attributes.add(classAtt);
        Instances res = new Instances("trainingData", attributes, capacity);
        res.setClass(classAtt);

        return res;

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return m_classifier.toString();
    }

}
