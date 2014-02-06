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

import java.util.Iterator;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class ConcatVectorDataList implements VectorDataList {

    private final VectorDataList m_vecList2;
    private final VectorDataList m_vecList1;

    public ConcatVectorDataList(VectorDataList vecList1, VectorDataList vecList2) {
        m_vecList1 = vecList1;
        m_vecList2 = vecList2;

    }

    /**
     * {@inheritDoc}
     */
    public double[] getVector(int idx) {
        if (idx >= m_vecList1.numVectors()) {
            return m_vecList2.getVector(idx - m_vecList1.numVectors());
        } else {
            return m_vecList1.getVector(idx);
        }
    }

    /**
     * {@inheritDoc}
     */
    public int getClusterIdx(int vecIdx) {
        if (vecIdx >= m_vecList1.numVectors()) {
            return m_vecList2.getClusterIdx(vecIdx - m_vecList1.numVectors());
        } else {
            return m_vecList1.getClusterIdx(vecIdx);
        }
    }

    /**
     * {@inheritDoc}
     */
    public void setClusterIdx(int vecIdx, int clusterIdx) {
        if (vecIdx >= m_vecList1.numVectors()) {
            m_vecList2.setClusterIdx(vecIdx - m_vecList1.numVectors(),
                    clusterIdx);
        } else {
            m_vecList1.setClusterIdx(vecIdx, clusterIdx);
        }

    }

    /**
     * {@inheritDoc}
     */
    public int numVectors() {
        return m_vecList1.numVectors() + m_vecList2.numVectors();
    }

    /**
     * {@inheritDoc}
     */
    public int numVectors(int clustIdx) {
        return m_vecList1.numVectors(clustIdx)
                + m_vecList2.numVectors(clustIdx);
    }

    /**
     * {@inheritDoc}
     */
    public int numClusters() {
        // TODO Auto-generated method stub
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    public int numFeatures() {
        return m_vecList1.numFeatures();
    }

    /**
     * {@inheritDoc}
     */
    public Iterator<Integer> iterator(int clustIdx) {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     * {@inheritDoc}
     */
    public Iterator<double[]> iterator() {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     * {@inheritDoc}
     */
    public void setWeight(int vecIdx, double weight) {
        // TODO Auto-generated method stub

    }

    /**
     * {@inheritDoc}
     */
    public double weight(int vecIdx) {
        // TODO Auto-generated method stub
        return 0;
    }

}
