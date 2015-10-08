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

import org.knime.knip.core.KNIPGateway;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.numeric.RealType;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public abstract class ContourDataExtractor extends AbstractVectorDataList {

	private ContourDataGrid m_cDataGrid;

	private int[] m_translations;

	private int[] m_permutation;

	protected int CENTER_COL;

	/**
	 * 
	 */
	public void extractContourData(ContourDataGrid cDataGrid) {
		m_cDataGrid = cDataGrid;
		CENTER_COL = cDataGrid.width() / 2;
		m_translations = new int[cDataGrid.totalLength()];
		m_permutation = new int[cDataGrid.totalLength()];
		for (int i = 0; i < m_permutation.length; i++) {
			m_permutation[i] = i;
		}
		extractContourData(m_translations, m_permutation);
	}

	/**
	 * @return
	 */
	protected ContourDataGrid contourDataGrid() {
		return m_cDataGrid;
	}

	/**
	 * @param cDataGrid
	 * @param translations
	 * @param permutation
	 */
	protected abstract void extractContourData(int[] translations, int[] permutation);

	/**
	 * {@inheritDoc}
	 */
	@Override
	public double[] getVector(int idx) {
		return m_cDataGrid.get(idx, CENTER_COL + m_translations[idx]);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int numVectors() {
		return m_cDataGrid.totalLength();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int numFeatures() {
		return m_cDataGrid.numFeatures();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Iterator<double[]> iterator() {
		List<double[]> tmp = new ArrayList<double[]>(numVectors());
		for (int i = 0; i < numVectors(); i++) {
			tmp.add(getVector(i));
		}
		return tmp.iterator();
	}

	/**
	 * @return the number of samples included
	 */
	public int numSamples() {
		return m_cDataGrid.numSamples();
	}

	/**
	 * @param sampleIdx
	 * @return the length of the sample
	 */
	public int getSampleLength(int sampleIdx) {
		return m_cDataGrid.getSampleLength(sampleIdx);
	}

	/**
	 * @param sampleIdx
	 * @return the accumulated length of all the samples to the specified sample
	 *         index (without the sample of <code>sampleIdx</code> itself!)
	 */
	public int getAccumulatedSampleLength(int sampleIdx) {
		return m_cDataGrid.getAccumulatedSampleLength(sampleIdx);
	}

	/**
	 * @param rowIdx
	 * @return the sample index to the global row index
	 */
	public int getSampleIndex(int rowIdx) {
		return m_cDataGrid.getSampleIndex(rowIdx);
	}

	public List<double[]> nonContourVectors() {
		List<double[]> res = new ArrayList<double[]>(m_cDataGrid.numVectors() - m_cDataGrid.totalLength());
		for (int i = 0; i < m_cDataGrid.width(); i++) {
			for (int j = 0; j < m_cDataGrid.totalLength(); j++) {
				if (i != CENTER_COL + m_translations[j]) {
					res.add(m_cDataGrid.get(j, i));
				}
			}
		}
		return res;
	}

	public double getCentrality() {
		double res = 0;
		for (int i = 0; i < m_translations.length; i++) {
			res += m_translations[i] * weight(i);
		}
		return res / m_translations.length;
	}

	public double getContinuity() {
		// TODO: ignore samples transitions and zero-weighted lines in
		// continuity determination
		double res = 0;
		for (int i = 1; i < m_translations.length; i++) {
			res = Math.max(res, Math.abs(m_translations[i] - m_translations[i - 1]) * weight(i) * weight(i - 1));
		}
		res = Math.max(res, Math.abs(m_translations[0] - m_translations[m_translations.length - 1]) * weight(0)
				* weight(m_translations.length - 1));
		return res;
	}

	public RandomAccessibleInterval<LabelingType<Integer>> clusterDistrLabeling(Integer bgCluster) {
		// read labeling mapping and create labeling
		long[] dims = new long[] { m_cDataGrid.width(), numVectors() };
		RandomAccessibleInterval<LabelingType<Integer>> res = (RandomAccessibleInterval<LabelingType<Integer>>) KNIPGateway
				.ops().create().imgLabeling(dims);

		RandomAccess<LabelingType<Integer>> ra = res.randomAccess();
		for (int h = 0; h < res.dimension(1); h++) {
			int clustIdx;
			Integer label;
			if ((clustIdx = getClusterIdx(h)) != bgCluster && weight(h) > 0) {
				label = clustIdx;
			} else {
				label = null;
			}
			ra.setPosition(h, 1);
			for (int w = 0; w < res.dimension(0); w++) {
				ra.setPosition(w, 0);
				if (label == null)
					ra.get().clear();
				else
					ra.get().add(label);
			}
		}
		return res;
	}

	/*
	 * Debugging
	 */
	public <T extends RealType<T>> Img<T> transformContourImage(Img<T> contourImg) {

		Img<T> res = contourImg.copy();
		RandomAccess<T> resRA = res.randomAccess();
		Cursor<T> srcCur = contourImg.localizingCursor();

		while (srcCur.hasNext()) {
			srcCur.fwd();
			resRA.setPosition((srcCur.getIntPosition(0) - m_translations[srcCur.getIntPosition(1)] + res.dimension(0))
					% res.dimension(0), 0);
			resRA.setPosition(m_permutation[srcCur.getIntPosition(1)], 1);
			resRA.get().set(srcCur.get());
		}

		return res;
	}

}
