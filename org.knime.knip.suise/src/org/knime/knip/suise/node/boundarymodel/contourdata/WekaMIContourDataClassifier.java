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
import weka.core.Attribute;
import weka.core.Capabilities.Capability;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

/**
 * 
 * Transforms the Contour data to a multi-instance data set and trains a mi
 * classifier. The contour data has class index 1!
 * 
 * @author hornm, University of Konstanz
 */
public class WekaMIContourDataClassifier implements ContourDataClassifier {

    /* multi-instance classifier */
    private final Classifier m_classifier;

    private Instances m_data;

    /**
     * @param classifier a multi-instance classifier
     */
    public WekaMIContourDataClassifier(Classifier classifier) {
        if (!classifier.getCapabilities()
                .handles(Capability.ONLY_MULTIINSTANCE)) {
            throw new IllegalArgumentException(
                    "Classifier can not handle multi-instances");
        }
        m_classifier = classifier;

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void buildClassifier(ContourDataGrid cData, VectorDataList bgData)
            throws Exception {

        // transform input data to weka mi-instances
        m_data =
                initDataset(cData.numFeatures(), 2, cData.totalLength()
                        + bgData.numVectors(), cData.width());

        for (int r = 0; r < cData.totalLength(); r++) {
            Instances bagData =
                    new Instances(m_data.attribute(1).relation(), cData.width());
            for (int c = 0; c < cData.width(); c++) {
                int vecIdx = cData.getVectorIdx(c, r);
                Instance inst =
                        new DenseInstance(cData.weight(vecIdx),
                                cData.getVector(vecIdx));
                inst.setDataset(bagData);
                bagData.add(inst);
            }
            int value = m_data.attribute(1).addRelation(bagData);
            Instance newBag = new DenseInstance(3);
            newBag.setValue(0, r); // bag id
            newBag.setValue(2, 1); // class attribute
            newBag.setValue(1, value);
            newBag.setWeight(1);
            newBag.setDataset(m_data);
            m_data.add(newBag);
        }

        for (int i = 0; i < bgData.numVectors(); i++) {
            Instances bagData =
                    new Instances(m_data.attribute(1).relation(), cData.width());
            Instance inst =
                    new DenseInstance(bgData.weight(i), bgData.getVector(i));
            inst.setDataset(bagData);
            bagData.add(inst);
            int value = m_data.attribute(1).addRelation(bagData);
            Instance newBag = new DenseInstance(3);
            newBag.setValue(0, cData.totalLength() + i);
            newBag.setValue(2, 0);
            newBag.setValue(1, value);
            newBag.setWeight(1);
            newBag.setDataset(m_data);
            m_data.add(newBag);
        }

        m_classifier.buildClassifier(m_data);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double contourProbability(double[] inst) throws Exception {
        Instances bagData = new Instances(m_data.attribute(1).relation(), 1);
        Instance i = new DenseInstance(1, inst);
        i.setDataset(bagData);

        bagData.add(i);

        Instance bag = new DenseInstance(3);
        bag.setDataset(m_data);
        int val = bag.attribute(1).addRelation(bagData);
        bag.setValue(1, val);

        return m_classifier.distributionForInstance(bag)[1];

    }

    /*
     * Helper method
     */
    private Instances initDataset(int numFeatures, int numClasses, int numBags,
            int capacityPerBag) {
        ArrayList<Attribute> attributes =
                new ArrayList<Attribute>(numFeatures + 5);
        for (int i = 0; i < numFeatures; i++) {
            attributes.add(new Attribute("att" + i));
        }
        Instances bagData = new Instances("bagData", attributes, numBags);

        ArrayList<String> classNames = new ArrayList<String>();
        for (int i = 0; i < numClasses; i++) {
            classNames.add("class" + i);
        }
        Attribute attClass = new Attribute("class", classNames);

        ArrayList<String> bagID = new ArrayList<String>(numBags);
        for (int i = 0; i < numBags; i++) {
            bagID.add("bag#" + i);
        }
        Attribute attBagIndex = new Attribute("bag index", bagID);
        ArrayList<Attribute> attInfo = new ArrayList<Attribute>(3);
        attInfo.add(attBagIndex);
        attInfo.add(new Attribute("bag", bagData)); // relation-valued attribute
        attInfo.add(attClass);
        Instances res = new Instances("Multi-Instance-Dataset", attInfo, 0);
        res.setClassIndex(res.numAttributes() - 1);

        return res;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] contourProbDistribution(double[] inst) throws Exception {
        return new double[]{contourProbability(inst)};
    }
}
