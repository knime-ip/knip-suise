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

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.AbstractClassifier;
import weka.classifiers.UpdateableClassifier;
import weka.core.Capabilities;
import weka.core.Capabilities.Capability;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

/**
 * 
 * @author hornm, University of Konstanz
 */
public class IntervalRule extends AbstractClassifier implements
        UpdateableClassifier {

    private Instance m_min;
    private Instance m_max;

    public IntervalRule() {
    }

    /**
     * Copy constructor.
     * 
     * @param ir
     */
    public IntervalRule(IntervalRule ir) {
        m_min = new DenseInstance(1.0, ir.m_min.toDoubleArray());
        m_max = new DenseInstance(1.0, ir.m_max.toDoubleArray());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void buildClassifier(Instances data) throws Exception {
        // can classifier handle the data?
        getCapabilities().testWithFail(data);

        for (Instance i : data) {
            updateClassifier(i);
        }

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void updateClassifier(Instance instance) throws Exception {
        if (m_min == null || m_max == null) {
            m_min = new DenseInstance(1.0, instance.toDoubleArray());
            m_max = new DenseInstance(1.0, instance.toDoubleArray());
        } else {
            for (int a = 0; a < instance.numAttributes() - 1; a++) {
                m_min.setValue(a, Math.min(instance.value(a), m_min.value(a)));
                m_max.setValue(a, Math.max(instance.value(a), m_max.value(a)));
            }
        }

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] distributionForInstance(Instance instance) throws Exception {
        for (int a = 0; a < instance.numAttributes() - 1; a++) {
            if (instance.value(a) < m_min.value(a)
                    || instance.value(a) > m_max.value(a)) {
                return new double[] { 1, 0 };
            }
        }
        return new double[] { 0, 1 };

        // int fit = 0;
        // for (int a = 0; a < instance.numAttributes() - 1; a++) {
        // if (instance.value(a) >= m_min.value(a)
        // || instance.value(a) <= m_max.value(a)) {
        // fit++;
        // }
        // }
        // double r = (double) fit / (instance.numAttributes() - 1);
        // return new double[] { 1 - r, r };

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Capabilities getCapabilities() {
        Capabilities result = super.getCapabilities();

        // attributes
        result.enable(Capability.NUMERIC_ATTRIBUTES);

        // class
        result.disableAllClasses();
        result.disableAllClassDependencies();
        result.enable(Capability.BINARY_CLASS);

        return result;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return "min=" + m_min + ";max=" + m_max;
    }

}
