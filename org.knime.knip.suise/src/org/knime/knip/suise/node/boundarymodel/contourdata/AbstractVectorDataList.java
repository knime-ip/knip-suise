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
import java.util.HashMap;
import java.util.Iterator;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public abstract class AbstractVectorDataList implements VectorDataList {

    private HashMap<Integer, Double> m_weightMap;
    private HashMap<Integer, Integer> m_vectorClusterMap;
    private ArrayList<ArrayList<Integer>> m_vecClusterLists;
    private int m_numCluster = -1;

    public AbstractVectorDataList() {
        m_weightMap = new HashMap<Integer, Double>();
        m_vectorClusterMap = new HashMap<Integer, Integer>();
        m_vecClusterLists = new ArrayList<ArrayList<Integer>>();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getClusterIdx(int vecIdx) {
        Integer res = m_vectorClusterMap.get(vecIdx);
        if (res == null) {
            return 0;
        } else {
            return res;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setClusterIdx(int vecIdx, int clusterIdx) {
        double[] vec = getVector(vecIdx);
        if (vec != null) {
            if (clusterIdx != 0) {
                m_vectorClusterMap.put(vecIdx, clusterIdx);
            } else {
                m_vectorClusterMap.remove(vecIdx);
            }
        }
        m_numCluster = -1;

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numVectors(int clustIdx) {
        finalizeClusterStructure();
        return m_vecClusterLists.get(clustIdx).size();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numClusters() {
        finalizeClusterStructure();
        return m_numCluster;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Iterator<Integer> iterator(int clustIdx) {
        finalizeClusterStructure();
        return m_vecClusterLists.get(clustIdx).iterator();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double weight(int vecIdx) {
        Double res = m_weightMap.get(vecIdx);
        if (res == null) {
            return 1.0;
        } else {
            return res;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setWeight(int vecIdx, double weight) {
        if (getVector(vecIdx) != null) {
            if (weight == 1.0) {
                if (m_weightMap.get(vecIdx) != null) {
                    m_weightMap.remove(vecIdx);
                }
            } else {
                m_weightMap.put(vecIdx, weight);
            }
        }

    }

    private void finalizeClusterStructure() {
        if (m_numCluster == -1) {
            m_vecClusterLists.clear();
            int clustIdx;
            for (int vecIdx = 0; vecIdx < numVectors(); vecIdx++) {
                clustIdx = getClusterIdx(vecIdx);
                m_numCluster = Math.max(clustIdx + 1, m_numCluster);
                while (m_vecClusterLists.size() < m_numCluster) {
                    m_vecClusterLists.add(new ArrayList<Integer>(numVectors()));
                }
                m_vecClusterLists.get(clustIdx).add(vecIdx);
            }
        }
        // TODO: trim to size
    }
}
