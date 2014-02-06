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
import java.util.List;

import weka.core.Utils;
import weka.core.matrix.Matrix;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class ContourCluster implements Iterable<int[]> {

    public enum ClusterScore {
        ESS, EUCL_DIST, LMDL;
    }

    private List<int[]> m_samples;
    private double m_score = Double.NaN;
    private double m_formerScore = Double.NaN;
    private ClusterScore m_formerType = ClusterScore.EUCL_DIST;
    private double[] m_centroid = null;
    private ContourDataGrid m_cdata;

    public ContourCluster(ContourDataGrid cdata) {
        m_cdata = cdata;
        m_samples = new ArrayList<int[]>(m_cdata.totalLength());
    }

    public ContourCluster(ContourDataGrid data, double[] centroid) {
        // TODO:
    }

    public void addSample(int colIdx, int rowIdx) {
        m_samples.add(new int[] { colIdx, rowIdx });
        m_formerScore = m_score;
        m_score = Double.NaN;
        m_centroid = null;
    }

    public void addCluster(ContourCluster c) {
        m_samples.addAll(c.m_samples);
        m_formerScore = Double.NaN;
        m_score = Double.NaN;
        m_centroid = null;
    }

    public void removeLastSample() {
        m_samples.remove(m_samples.size() - 1);
        m_score = m_formerScore;
        m_formerScore = Double.NaN;
        m_centroid = null;
    }

    public double getClusterScore(ClusterScore cs) {
        if (Double.isNaN(m_score) || !m_formerType.equals(cs)) {
            m_score = calcClusterScore(cs);
        }
        return m_score;
    }

    public int size() {
        return m_samples.size();
    }

    public int[] getSample(int index) {
        return m_samples.get(index);
    }

    public void clear() {
        m_samples.clear();
        m_centroid = null;
        m_score = Double.NaN;
        m_formerScore = Double.NaN;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Iterator<int[]> iterator() {
        return m_samples.iterator();
    }

    private double calcClusterScore(ClusterScore cs) {

        switch (cs) {
        case EUCL_DIST:
            // distance to the first cluster element as score
            return dist(m_cdata.get(m_samples.get(0)[0], m_samples.get(0)[1]),
                    m_cdata.get(m_samples.get(m_samples.size() - 1)[0],
                            m_samples.get(m_samples.size() - 1)[1]));
        case ESS:
            // error sum-of-squares as score
            return calcESS();
        case LMDL:
            // lossy minimum description length
            return calcLMDL();

        }

        return 0;

    }

    public double clusterDist(ContourCluster c2) {
        return dist(centroid(), c2.centroid());
    }

    /*
     * Euclidean distance
     */
    private double dist(double[] v1, double[] v2) {
        double res = 0;
        for (int i = 0; i < v2.length; i++) {
            res += Math.pow(v1[i] - v2[i], 2);
        }
        return Math.sqrt(res);
    }

    private double[] centroid() {

        if (m_centroid == null) {

            m_centroid = new double[m_cdata.numFeatures()];
            for (int i = 0; i < m_samples.size(); i++) {
                double[] vec = m_cdata.get(m_samples.get(i)[0],
                        m_samples.get(i)[1]);
                for (int j = 0; j < m_cdata.numFeatures(); j++) {
                    m_centroid[j] += vec[j];
                }
            }

            for (int j = 0; j < m_cdata.numFeatures(); j++) {
                m_centroid[j] /= m_samples.size();
            }
        }

        return m_centroid;

    }

    public double calcCentroidDist() {
        // TODO:
        centroid();
        return 0;
    }

    /** calculated error sum-of-squares for instances wrt centroid **/
    private double calcESS() {

        double[] centroid = centroid();

        double fESS = 0;
        double[] vec;
        for (int i = 0; i < m_samples.size(); i++) {
            vec = m_cdata.get(m_samples.get(i)[0], m_samples.get(i)[1]);
            fESS += dist(centroid, vec);
        }
        return fESS / m_samples.size();
    } // calcESS

    private double calcLMDL() {

        // create sample matrix
        Matrix V = new Matrix(m_cdata.numFeatures(), m_samples.size());

        double[] vec;
        for (int i = 0; i < m_samples.size(); i++) {
            vec = m_cdata.get(m_samples.get(i)[0], m_samples.get(i)[1]);
            for (int j = 0; j < vec.length; j++) {
                V.set(j, i, vec[j]);
            }
        }

        double epsilon = 5;

        // estimate of the covariance matrix
        Matrix W = V.times(V.transpose());

        W.times(m_cdata.numFeatures() / (epsilon * epsilon * m_samples.size()));

        W = Matrix.identity(m_cdata.numFeatures(), m_cdata.numFeatures()).plus(
                W);

        return Utils.log2(W.det()) * (m_cdata.numFeatures() + m_samples.size())
                / 2;

    }

}
