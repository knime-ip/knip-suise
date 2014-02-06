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
import java.util.Collections;
import java.util.HashMap;

/**
 * Permutohedral Lattice for High-Dimensional Filtering
 * 
 * @see @ARTICLE{Adams2010, author = {Adams, Andrew and Baek, Jongmin and Davis,
 *      Myers Abraham}, title = {Fast High-Dimensional Filtering Using the
 *      Permutohedral Lattice}, journal = {Computer Graphics Forum}, year =
 *      {2010}, volume = {29}, pages = {753--762}, number = {2}, doi =
 *      {10.1111/j.1467-8659.2009.01645.x}, file =
 *      {Adams2010.pdf:Adams2010.pdf:PDF}, issn = {1467-8659}, keywords = {I.4.3
 *      [Image Processing and Computer Vision]: Enhancement?Filtering},
 *      publisher = {Blackwell Publishing Ltd}, url =
 *      {http://dx.doi.org/10.1111/j.1467-8659.2009.01645.x} }
 * 
 *      Implementation according to the c-implementation of Adams2010
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class PermutohedralLattice {

    private short[] m_canonical;

    private float[] m_scaleFactor;

    private short[] m_greedy;

    private byte[] m_rank;

    private float[] m_barycentric;

    private int m_nReplay;

    private short[] m_key;

    private float[] m_elevated;

    private int m_d;

    private int m_vd;

    private HashMap<Integer, Integer> m_keyIndexTable;

    private short[][] m_keys;

    private float[][] m_values;

    private float[][] m_tmpValues;

    // weights for each simplex vertex to be replayed for the slicing
    private float[] m_replay;

    /**
     * @param d dimensionality of the key vectors
     * @param vd dimensionality of the value vectors
     * @param nData number of points in the input
     */
    public PermutohedralLattice(int d, int vd, int nData) {

        m_elevated = new float[d + 1];
        m_scaleFactor = new float[d];

        // coordinates of the canonical simplex
        m_canonical = new short[(d + 1) * (d + 1)];

        m_greedy = new short[d + 1];
        m_rank = new byte[d + 1];
        m_barycentric = new float[d + 2];
        m_replay = new float[nData * (d + 1)];
        m_nReplay = 0;
        m_key = new short[d + 1];

        m_d = d;
        m_vd = vd;

        m_keyIndexTable = new HashMap<Integer, Integer>(nData * (d + 1));
        m_values = new float[nData * (d + 1)][];
        m_tmpValues = new float[(nData * (d + 1))][];
        m_keys = new short[nData * (d + 1)][];

        // compute the coordinates of the canonical simplex, in which the
        // difference between a contained point and the zero remainder vertex is
        // always in ascending order
        for (int i = 0; i <= d; i++) {
            for (int j = 0; j <= d - i; j++)
                m_canonical[i * (d + 1) + j] = (short)i;
            for (int j = d - i + 1; j <= d; j++)
                m_canonical[i * (d + 1) + j] = (short)(i - (d + 1));
        }

        // compute parts of the rotation matrix E (orthogonal basis)
        for (int i = 0; i < d; i++) {
            // the diagonal entries for normalization
            m_scaleFactor[i] = (float)(1.0 / Math.sqrt((i + 1) * (i + 2)));

            /*
             * We presume that the user would like to do a Gaussian blur of
             * standard deviation 1 in each dimension (or a total variance of d,
             * summed over dimensions.) Because the total variance of the blur
             * performed by this algorithm is not d, we must scale the space to
             * offset this.
             * 
             * The total variance of the algorithm is (See pg.6 and 10 of
             * paper): [variance of splatting] + [variance of blurring] +
             * [variance of splatting] = d(d+1)(d+1)/12 + d(d+1)(d+1)/2 +
             * d(d+1)(d+1)/12 = 2d(d+1)(d+1)/3.
             * 
             * So we need to scale the space by (d+1)sqrt(2/3).
             */
            m_scaleFactor[i] *= (d + 1) * Math.sqrt(2.0 / 3);
        }

    }

    /**
     * Performs splatting with given position and value vectors
     * 
     * @param position the position with d entries
     * @param value the value with vd entries
     */
    public void splat(double[] position, double[] value) {
        // first rotate position into the (d+1)-dimensional hyperplane
        m_elevated[m_d] =
                (float)(-m_d * position[m_d - 1] * m_scaleFactor[m_d - 1]);
        for (int i = m_d - 1; i > 0; i--) {
            m_elevated[i] =
                    (float)(m_elevated[i + 1] - i * position[i - 1]
                            * m_scaleFactor[i - 1] + (i + 2) * position[i]
                            * m_scaleFactor[i]);
        }
        m_elevated[0] =
                m_elevated[1] + (float)(2 * position[0] * m_scaleFactor[0]);

        // prepare to find the closest lattice points
        float scale = 1.0f / (m_d + 1);
        Arrays.fill(m_rank, (byte)0);

        // greedily search for the closest zero-colored lattice point
        int sum = 0;
        for (int i = 0; i <= m_d; i++) {
            float v = m_elevated[i] * scale;
            double up = Math.ceil(v) * (m_d + 1);
            double down = Math.floor(v) * (m_d + 1);

            if (up - m_elevated[i] < m_elevated[i] - down)
                m_greedy[i] = (short)up;
            else
                m_greedy[i] = (short)down;

            sum += m_greedy[i];
        }
        sum /= m_d + 1;

        // rank differential to find the permutation between this simplex and
        // the canonical one.
        // (See pg. 3-4 in paper.)
        // memset(myrank, 0, sizeof(char)*(d+1));
        for (int i = 0; i < m_d; i++)
            for (int j = i + 1; j <= m_d; j++)
                if (m_elevated[i] - m_greedy[i] < m_elevated[j] - m_greedy[j])
                    m_rank[i]++;
                else
                    m_rank[j]++;

        if (sum > 0) {
            // sum too large - the point is off the hyperplane.
            // need to bring down the ones with the smallest differential
            for (int i = 0; i <= m_d; i++) {
                if (m_rank[i] >= m_d + 1 - sum) {
                    m_greedy[i] -= m_d + 1;
                    m_rank[i] += sum - (m_d + 1);
                } else
                    m_rank[i] += sum;
            }
        } else if (sum < 0) {
            // sum too small - the point is off the hyperplane
            // need to bring up the ones with largest differential
            for (int i = 0; i <= m_d; i++) {
                if (m_rank[i] < -sum) {
                    m_greedy[i] += m_d + 1;
                    m_rank[i] += (m_d + 1) + sum;
                } else
                    m_rank[i] += sum;
            }
        }

        // Compute barycentric coordinates (See pg.10 of paper.)
        Arrays.fill(m_barycentric, 0);
        // memset(barycentric, 0, sizeof(float)*(d+2));
        for (int i = 0; i <= m_d; i++) {
            m_barycentric[m_d - m_rank[i]] +=
                    (m_elevated[i] - m_greedy[i]) * scale;
            m_barycentric[m_d + 1 - m_rank[i]] -=
                    (m_elevated[i] - m_greedy[i]) * scale;
        }
        m_barycentric[0] += 1.0f + m_barycentric[m_d + 1];

        // Splat the value into each vertex of the simplex, with barycentric
        // weights.
        for (int remainder = 0; remainder <= m_d; remainder++) {
            // Compute the location of the lattice point explicitly (all but the
            // last coordinate - it's redundant because they sum to zero)
            for (int i = 0; i < m_d; i++) {
                m_key[i] =
                        (short)(m_greedy[i] + m_canonical[remainder * (m_d + 1)
                                + m_rank[i]]);
            }

            // Retrieve pointer to the value at this vertex.
            // float * val = hashTable.lookup(m_key, true);

            // Accumulate values with barycentric weight.

            int hash = hashCode(m_key);
            Integer index = m_keyIndexTable.get(hash);

            float[] val;
            float[] tmpVal;
            if (index == null) {
                val = new float[m_vd];
                tmpVal = new float[m_vd];
                m_keyIndexTable.put(hash, m_nReplay);
                m_keys[m_nReplay] = m_key.clone();
            } else {
                val = m_values[index];
                tmpVal = m_tmpValues[index];
                m_keys[m_nReplay] = m_keys[index];
            }
            m_values[m_nReplay] = val;
            m_tmpValues[m_nReplay] = tmpVal;

            for (int i = 0; i < m_vd; i++)
                val[i] += m_barycentric[remainder] * value[i];

            // Record this interaction to use later when slicing
            // re.offset = val - hashTable.getValues();

            m_replay[m_nReplay] = m_barycentric[remainder];
            m_nReplay++;

        }

    }

    /**
     * 
     */
    public void beginSplatValue() {
        m_nReplay = 0;
        for (int i = 0; i < m_values.length; i++) {
            Arrays.fill(m_values[i], 0);
            Arrays.fill(m_tmpValues[i], 0);
        }
    }

    /**
     * "Re-splates" the values using the already calculated barycentric weights
     * (splat(double[], double[]) has to be called once, first)
     * 
     * @param value the value(s)
     */
    public void splatValue(double[] value) {
        for (int i = 0; i <= m_d; i++) {
            float[] values = m_values[m_nReplay];
            for (int j = 0; j < m_vd; j++) {
                // col[j] += r.weight * base[r.offset + j];
                values[j] += (float)(m_replay[m_nReplay] * value[j]);
                // col[j] = weights[i] * values[j];
            }
            m_nReplay++;
        }
    }

    /**
     * Prepare for slicing
     */
    public void beginSlice() {
        m_nReplay = 0;
    }

    /**
     * Performs slicing out of position vectors. Note that the barycentric
     * weights and the simplex containing each position vector were calculated
     * and stored in the splatting step. We may reuse this to accelerate the
     * algorithm. (See pg. 6 in paper.)
     */
    public void slice(double[] col) {
        // float *base = hashTable.getValues();
        Arrays.fill(col, 0);
        for (int i = 0; i <= m_d; i++) {
            float[] values = m_values[m_nReplay];
            for (int j = 0; j < m_vd; j++) {
                // col[j] += r.weight * base[r.offset + j];
                col[j] += m_replay[m_nReplay] * values[j];
            }
            m_nReplay++;
        }

    }

    /** Performs a Gaussian blur along each projected axis in the hyperplane. */
    public void blur() {

        // Prepare arrays
        short[] neighbor1 = new short[m_d + 1];
        short[] neighbor2 = new short[m_d + 1];
        // float[] newValue = new float[m_vd*hashTable.size()];
        float[] oldValue;
        // float [] hashTableBase = oldValue;

        float[] zero = new float[m_vd];
        // for (int k = 0; k < m_vd; k++) zero[k] = 0;

        // For each of d+1 axes,
        for (int j = 0; j <= m_d; j++) {
            // printf(" %d", j);fflush(stdout);

            ArrayList<Integer> tmpList =
                    new ArrayList<Integer>(m_keyIndexTable.values());
            Collections.sort(tmpList);

            // For each vertex in the lattice,
            for (Integer index : tmpList) { // blur
                                            // point
                                            // i
                                            // in
                                            // dimension
                                            // j
                short[] key = m_keys[index]; // keys to current
                                             // vertex

                // short[] key = hashTable.getKeys() + i*(m_d);

                for (int k = 0; k < m_d; k++) {
                    neighbor1[k] = (short)(key[k] + 1);
                    neighbor2[k] = (short)(key[k] - 1);
                }
                neighbor1[j] = (short)(key[j] - m_d);
                neighbor2[j] = (short)(key[j] + m_d); // keys to the neighbors
                                                      // along the given axis.

                oldValue = m_values[index];

                // float *oldVal = oldValue + i*m_vd;
                // float *newVal = newValue + i*m_vd;

                // float *vm1, *vp1;

                float[] vm1, vp1;

                // vm1 = hashTable.lookup(neighbor1, false);
                Integer neighborIndex =
                        m_keyIndexTable.get(hashCode(neighbor1));// look
                // up
                // first
                // neighbor

                // if (vm1!=null) vm1 = vm1 - hashTableBase + oldValue;
                // else vm1 = zero;
                if (neighborIndex == null) {
                    vm1 = zero;
                } else {
                    vm1 = m_values[neighborIndex];
                }

                // vp1 = hashTable.lookup(neighbor2, false);
                neighborIndex = m_keyIndexTable.get(hashCode(neighbor2));// look
                                                                         // up
                                                                         // second
                // neighbor
                // if (vp1!=null) vp1 = vp1 - hashTableBase + oldValue;
                // else vp1 = zero;
                if (neighborIndex == null) {
                    vp1 = zero;
                } else {
                    vp1 = m_values[neighborIndex];
                }

                // Mix values of the three vertices
                float[] newValue = m_tmpValues[index];
                for (int k = 0; k < m_vd; k++)
                    newValue[k] =
                            (0.25f * vm1[k] + 0.5f * oldValue[k] + 0.25f * vp1[k]);

            }

            float[][] tmp = m_tmpValues;
            m_tmpValues = m_values;
            m_values = tmp;

        }

        // the freshest data is now in oldValue, and newValue is ready to be
        // written over

        // depending where we ended up, we may have to copy data
        // if (oldValue != hashTableBase) {
        // memcpy(hashTableBase, oldValue, hashTable.size()*vd*sizeof(float));
        // delete oldValue;
        // } else {
        // delete newValue;
        // }
        // printf("\n");
        //
        // delete zero;
        // delete neighbor1;
        // delete neighbor2;
    }

    public int d() {
        return m_d;
    }

    private int hashCode(short[] a) {
        if (a == null)
            return 0;

        int result = 1;
        for (int i = 0; i < a.length - 1; i++) {
            result = 31 * result + a[i];
        }
        return result;
    }

    public static void main(String[] args) {
        PermutohedralLattice lattice = new PermutohedralLattice(2, 1, 100);
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                lattice.splat(new double[]{j, i}, new double[]{1});
            }
        }
        lattice.blur();
        lattice.beginSlice();
        double[] val = new double[1];
        for (int i = 0; i < 100; i++) {
            lattice.slice(val);
            System.out.println(Arrays.toString(val));
        }
    }

}
