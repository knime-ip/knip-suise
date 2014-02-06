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
import java.util.Iterator;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class VectorDataListImpl extends AbstractVectorDataList {

    private ArrayList<double[]> m_vecList;

    /**
     * @param numClusters
     */
    public VectorDataListImpl() {
        m_vecList = new ArrayList<double[]>();
    }

    /**
     * Adds a vector clusterIdx 0 and weight 1.
     * 
     * @param vec
     */
    public void addVector(double[] vec) {
        addVector(vec, 0, 1.0);
    }

    /**
     * @param vec
     * @param clusterIdx
     */
    public void addVector(double[] vec, int clusterIdx) {
        addVector(vec, clusterIdx, 1.0);

    }

    /**
     * @param vec
     * @param clusterIdx
     * @param weight
     */
    public void addVector(double[] vec, int clusterIdx, double weight) {
        m_vecList.add(vec);
        setWeight(m_vecList.size() - 1, weight);
        setClusterIdx(m_vecList.size() - 1, clusterIdx);
    }

    /**
     * Trims the capacity.
     */
    public void trim() {
        m_vecList.trimToSize();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Iterator<double[]> iterator() {
        return m_vecList.iterator();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numVectors() {
        return m_vecList.size();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numFeatures() {
        return m_vecList.get(0).length;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getVector(int idx) {
        return m_vecList.get(idx);
    }

    public static void main(String[] args) {
        VectorDataListImpl vecList = new VectorDataListImpl();
        for (int i = 0; i < 1000; i++) {
            vecList.addVector(new double[10]);
        }

        for (int i = 0; i < 1000; i++) {
            vecList.addVector(new double[10], 4);
        }

        for (int i = 0; i < 1000; i++) {
            vecList.addVector(new double[10], 2, .5);
        }

        System.out.println("total num cluster:" + vecList.numClusters());
        for (int i = 0; i < vecList.numClusters(); i++) {
            System.out.println("vector in cluster " + i + ": "
                    + vecList.numVectors(i));

        }

        for (int i = 0; i < vecList.numVectors(); i += 500) {
            System.out.println("cluster and weight of vector " + i + ": "
                    + vecList.getClusterIdx(i) + " " + vecList.weight(i));
        }

        System.out.println("test finished");

    }

}
