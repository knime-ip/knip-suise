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
package org.knime.knip.suise.node.levelcca;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.converter.Converter;
import net.imglib2.converter.Converters;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.histogram.Integer1dBinMapper;
import net.imglib2.img.Img;
import net.imglib2.img.ImgView;
import net.imglib2.labeling.Labeling;
import net.imglib2.ops.operation.UnaryOperation;
import net.imglib2.ops.operation.iterableinterval.unary.MakeHistogram;
import net.imglib2.ops.operation.randomaccessibleinterval.unary.regiongrowing.CCA;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.LongType;

/**
 * Thresholds the input image at those intensity levels, which result in
 * different binary images and performs a connected component analysis on each.
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class LevelCCA<T extends IntegerType<T> & NativeType<T>> implements
		UnaryOperation<Img<T>, RandomAccessibleInterval<LabelingType<Integer>>> {

	private final long[][] m_structuringElement;

	private final boolean m_whiteBackground;

	private final int m_lowerBound;

	private final int m_upperBound;

	private final int m_stepSize;

	public LevelCCA(long[][] structuringElement, boolean whiteBackground) {
		this(structuringElement, whiteBackground, Integer.MIN_VALUE,
				Integer.MAX_VALUE, 1);

	}

	public LevelCCA(long[][] structuringElement, boolean whiteBackground,
			int lowerBound, int upperBound, int stepSize) {
		m_structuringElement = structuringElement;
		m_whiteBackground = whiteBackground;
		m_lowerBound = lowerBound;
		m_upperBound = upperBound;
		m_stepSize = stepSize;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public RandomAccessibleInterval<LabelingType<Integer>> compute(
			Img<T> input, RandomAccessibleInterval<LabelingType<Integer>> output) {
		CCA<BitType> cca = new CCA<BitType>(m_structuringElement, new BitType(
				m_whiteBackground));

		Histogram1d<T> hist = new Histogram1d<T>(new Integer1dBinMapper<T>(0,
				255, false));
		RandomAccess<LongType> raHist = hist.randomAccess();
		new MakeHistogram<T>((int) hist.getBinCount()).compute(input, hist);
		ThresholdConverter<T> c = new ThresholdConverter<T>(0);
		T type = input.firstElement().createVariable();
		long lastCountSum = 0;
		for (int i = 0; i < hist.getBinCount(); i += m_stepSize) {
			raHist.setPosition(i, 0);
			lastCountSum += raHist.get().getIntegerLong();
			if (i < m_lowerBound || i > m_upperBound || lastCountSum == 0) {
				continue;
			}
			lastCountSum = 0;
			hist.getCenterValue(i, type);
			c.setThreshold(type.getRealDouble());
			Img<BitType> tmp = null;
			try {
				tmp = new ImgView<BitType>(Converters.convert(
						(RandomAccessibleInterval<T>) input, c, new BitType()),
						input.factory().imgFactory(new BitType()));
			} catch (IncompatibleTypeException e) {
				throw new RuntimeException(e);
			}
			cca.compute(tmp, output);
		}

		return output;
	}

	/**
	 * {@inheritDoc}
	 */
	public UnaryOperation<Img<T>, RandomAccessibleInterval<LabelingType<Integer>>> copy() {
		return null;
	}

	private class ThresholdConverter<T extends RealType<T>> implements
			Converter<T, BitType> {

		private double m_threshold;

		private BitType m_type;

		public ThresholdConverter(double threshold) {
			m_threshold = threshold;
			m_type = new BitType();

		}

		public void setThreshold(double t) {
			m_threshold = t;
		}

		/**
		 * {@inheritDoc}
		 */
		public void convert(T input, BitType output) {
			output.set(input.getRealDouble() >= m_threshold);
		}
	}

}
