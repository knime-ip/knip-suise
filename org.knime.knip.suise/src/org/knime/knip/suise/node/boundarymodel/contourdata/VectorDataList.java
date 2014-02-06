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
 * Object representing a list of vectors, each assigned to a cluster.
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public interface VectorDataList extends Iterable<double[]> {

    /**
     * @param idx
     * @return the vector at the specified index
     */
    double[] getVector(int idx);

    /**
     * @param vecIdx
     * @return the cluster assignment of the vector at the specified indexF
     */
    int getClusterIdx(int vecIdx);

    /**
     * @param vecIdx
     * @param clusterIdx
     */
    void setClusterIdx(int vecIdx, int clusterIdx);

    /**
     * @return the number of total vectors
     */
    int numVectors();

    /**
     * 
     * @param clustIdx
     * @return number of vectors of the given cluster
     */
    int numVectors(int clustIdx);

    /**
     * @return number of clusters
     */
    int numClusters();

    /**
     * @return the length of the vectors
     */
    int numFeatures();

    /**
     * @param clustIdx
     * @return Iterator to iterate through the vector indices assigned to the
     *         specified cluster.
     */
    Iterator<Integer> iterator(int clustIdx);

    /**
     * {@inheritDoc} Note: iterates through ALL vectors.(also the disabled ones)
     */
    @Override
    public Iterator<double[]> iterator();

    /**
     * Sets the weight of the specified vector.
     * 
     * @param vecIdx
     * @param weight
     */
    void setWeight(int vecIdx, double weight);

    /**
     * @param vecIdx
     * @return the weight of the vector
     */
    double weight(int vecIdx);

}
