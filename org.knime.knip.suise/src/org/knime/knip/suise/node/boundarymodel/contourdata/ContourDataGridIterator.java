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
 * Iterates through the contour data, possibly consisting of multiple samples
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
class ContourDataGridIterator implements Iterator<double[]> {

    private double[][][] m_currentSample;
    private ContourDataGrid m_cdata;

    private int m_samplesIndex;
    private int m_NIdx; // length
    private int m_RIdx; // radius

    private int m_N; // length to the current sample
    private int m_globIdx;

    /**
     * @param contourDataGrid
     * 
     */
    public ContourDataGridIterator(ContourDataGrid contourDataGrid) {
        m_cdata = contourDataGrid;
        reset();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasNext() {
        return m_globIdx < m_cdata.width() * m_cdata.totalLength();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] next() {

        m_RIdx++;
        if (m_RIdx % m_cdata.width() == 0) {
            m_NIdx++;
            m_RIdx = 0;
        }
        if (m_NIdx - m_N >= m_currentSample.length) {
            m_N += m_currentSample.length;
            m_samplesIndex++;
            m_currentSample = m_cdata.getSample(m_samplesIndex);
        }
        m_globIdx++;
        return m_currentSample[m_NIdx - m_N][m_RIdx];

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException(
                "\"rmove() is not supported for the contour data iterator\"");

    }

    public int getPosX() {
        return m_RIdx;
    }

    public int getPosY() {
        return m_NIdx;
    }

    public int getCellSampleIndex() {
        return m_samplesIndex;
    }

    public void reset() {
        m_RIdx = -1;
        m_NIdx = -1;
        m_N = 0;
        m_samplesIndex = 0;
        m_currentSample = m_cdata.getSample(0);
        m_globIdx = 0;
    }

}
