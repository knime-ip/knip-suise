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
import java.util.List;
import java.util.Map;
import java.util.Vector;

import net.imglib2.img.Img;

import org.knime.knip.core.util.ShowInSameFrame.ImagePlaneProducer;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class ContourDataGrid extends AbstractVectorDataList {

    // list of samples with rows x columns x features (radius x length x
    // features)
    private List<double[][][]> m_cdata;

    private Map<Integer, Integer> m_rowToSampleMap;

    private List<Integer> m_accumulatedLengths;

    private int m_width;

    private int m_totalLength;

    private int m_numFeatures;

    // debugging
    protected Img m_contourDebugImage = null;

    protected ImagePlaneProducer m_imgProd = null;

    public ContourDataGrid() {
        m_cdata = new ArrayList<double[][][]>();

        m_width = 0;
        m_totalLength = 0;
        m_numFeatures = 0;

        m_rowToSampleMap = new HashMap<Integer, Integer>();
        m_accumulatedLengths = new Vector<Integer>();

    }

    /**
     * @param s NxRxA
     */
    public void addContourSample(double[][][] s) {
        m_accumulatedLengths.add(m_totalLength);
        for (int j = 0; j < s.length; j++) {
            m_rowToSampleMap.put(m_totalLength + j, m_cdata.size());
        }
        m_cdata.add(s);
        m_totalLength += s.length;
        m_width = s[0].length;
        m_numFeatures = s[0][0].length;

    }

    /**
     * 
     * @return The width of the contour data.
     */
    public int width() {
        return m_width;
    }

    /**
     * 
     * @return The global length of the contour data (over all samples).
     */
    public int totalLength() {
        return m_totalLength;
    }

    /**
     * @return the number of samples included
     */
    public int numSamples() {
        return m_cdata.size();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numFeatures() {
        return m_numFeatures;
    }

    /**
     * @param sampleIdx
     * @return the length of the sample
     */
    public int getSampleLength(int sampleIdx) {
        return m_cdata.get(sampleIdx).length;
    }

    /**
     * @param sampleIdx
     * @return the accumulated length of all the samples to the specified sample
     *         index (without the sample of <code>sampleIdx</code> itself!)
     */
    public int getAccumulatedSampleLength(int sampleIdx) {
        if (sampleIdx == m_accumulatedLengths.size()) {
            return totalLength();
        } else {
            return m_accumulatedLengths.get(sampleIdx);
        }
    }

    /**
     * @param rowIdx
     * @return the sample index to the global row index
     */
    public int getSampleIndex(int rowIdx) {
        return m_rowToSampleMap.get(rowIdx);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ContourDataGridIterator iterator() {
        return new ContourDataGridIterator(this);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getVector(int idx) {
        int row = idx / width();
        int col = idx % width();
        return get(row, col);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numVectors() {
        return width() * totalLength();
    }

    /**
     * For debugging purposes.
     * 
     * @param imgProd
     * @param contourDebugImage
     */
    public void setContourDebugImage(ImagePlaneProducer imgProd,
            Img contourDebugImage) {
        m_imgProd = imgProd;
        m_imgProd.updateConsumers(m_contourDebugImage);

    }

    /**
     * Calculates the vector index retrieved from the column and row index.
     * 
     * @param colIdx
     * @param rowIdx
     * @return
     */
    public int getVectorIdx(int colIdx, int rowIdx) {
        return rowIdx * width() + colIdx;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double weight(int vecIdx) {
        int col = vecIdx % width();
        // linear weight
        return (double)(4 * col * (width() - col))
                / (double)(width() * width());
    }

    // //////// Helper to access the contour data grid ///////////////

    /**
     * @return the feature vector at the global contour data position
     */
    protected double[] get(int rowIdx, int colIdx) {
        int sampleId = m_rowToSampleMap.get(rowIdx);
        return m_cdata.get(sampleId)[rowIdx
                - m_accumulatedLengths.get(sampleId)][colIdx];
    }

    /**
     * @return the feature vectors in the specified row
     */
    protected double[][] getRow(int rowIdx) {
        int sampleId = m_rowToSampleMap.get(rowIdx);
        return m_cdata.get(sampleId)[rowIdx
                - m_accumulatedLengths.get(sampleId)];
    }

    /**
     * @return the feature vectors of the sample
     */
    protected double[][][] getSample(int sampleIdx) {
        return m_cdata.get(sampleIdx);
    }

}
