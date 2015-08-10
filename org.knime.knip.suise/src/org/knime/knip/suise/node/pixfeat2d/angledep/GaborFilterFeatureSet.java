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
package org.knime.knip.suise.node.pixfeat2d.angledep;

import java.util.concurrent.ExecutorService;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.ops.img.UnaryConstantRightAssignment;
import net.imglib2.ops.img.UnaryOperationAssignment;
import net.imglib2.ops.operation.Operations;
import net.imglib2.ops.operation.complex.real.unary.ComplexImaginaryToRealAdapter;
import net.imglib2.ops.operation.complex.real.unary.ComplexRealToRealAdapter;
import net.imglib2.ops.operation.img.unary.ImgConvert;
import net.imglib2.ops.operation.iterable.unary.Sum;
import net.imglib2.ops.operation.iterableinterval.unary.MinMax;
import net.imglib2.ops.operation.real.binary.RealAdd;
import net.imglib2.ops.operation.real.unary.Normalize;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexDoubleType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.view.Views;

import org.knime.core.node.KNIMEConstants;
import org.knime.knip.base.KNIPConstants;
import org.knime.knip.core.KNIPGateway;
import org.knime.knip.core.ThreadPoolExecutorService;
import org.knime.knip.core.algorithm.convolvers.Convolver;
import org.knime.knip.core.algorithm.convolvers.ImgLib2FourierConvolver;
import org.knime.knip.core.algorithm.convolvers.filter.linear.Gabor;
import org.knime.knip.suise.data.feat.Dart;
import org.knime.knip.suise.data.feat.FeatImgGenerator;
import org.knime.knip.suise.data.feat.PixFeatureSet;

/**
 * 
 * @param <T>
 *            image type
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class GaborFilterFeatureSet<T extends RealType<T>, F extends RealType<F>>
		implements PixFeatureSet<T>, FeatImgGenerator<F> {

	private final ExecutorService m_executorService = new ThreadPoolExecutorService(
			KNIMEConstants.GLOBAL_THREAD_POOL.createSubPool(KNIPConstants.THREADS_PER_NODE));

	private Img<DoubleType>[] m_filters;

	private Cursor<DoubleType>[] m_filterCursors;

	private String[] m_names;

	private Img<F>[] m_convolvedImgs;

	private RandomAccess<F>[] m_convolvedImgRandAccess;

	private int m_currentAngID;

	private int m_numFeatures;

	private int m_halfNumAng;

	/* settings fields, to be serialized to restore this factory */
	private boolean m_precalcFeatures = true;

	private double[] m_scales;

	private double[] m_frequencies;

	private double[] m_elongations;

	private int m_radius;

	private int m_numAng;

	private long[] m_pos;

	private Img<T> m_img;

	private RandomAccess<T> m_tmpSrcRandAccess;

	private RandomAccessible<T> m_extendedImg;

	private F m_featureType;

	/**
	 * The user must not call this constructor.
	 */
	public GaborFilterFeatureSet() {
		// no-op
	}

	public GaborFilterFeatureSet(final double[] scales, final double[] frequencies, final double[] elongations,
			final int radius, final int numAng, F type) {
		this(scales, frequencies, elongations, radius, numAng, true, type);
	}

	/**
	 * @param scales
	 * @param frequencies
	 * @param elongations
	 * @param radius
	 * @param numAng
	 * @param precalcFeatures
	 * @param type
	 *            the feature type
	 * @param averageAngles
	 *            averages over the given number of angles, thus, updating the
	 *            angle has no effect
	 */
	public GaborFilterFeatureSet(final double[] scales, final double[] frequencies, final double[] elongations,
			final int radius, final int numAng, boolean precalcFeatures, F type) {
		super();

		m_scales = scales;
		m_frequencies = frequencies;
		m_elongations = elongations;
		m_numAng = numAng;
		m_precalcFeatures = precalcFeatures;
		m_radius = radius;
		m_featureType = type;
		init();
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	private void init() {

		m_halfNumAng = m_numAng / 2;

		// the total number of filters including the varied angles
		int numFilters = m_halfNumAng * m_scales.length * m_frequencies.length * m_elongations.length * 2;

		// create filter for the individual parameters
		m_filters = new Img[numFilters];
		m_names = new String[numFilters / m_halfNumAng];

		// TODO: ComplexImaginaryToRealAdapter CleanUp
		// calculating the filter for each angle, scale, frequency,
		// elongation
		// and symmetry (even or odd).
		UnaryOperationAssignment<ComplexDoubleType, DoubleType> realAssignment = new UnaryOperationAssignment(
				new ComplexRealToRealAdapter());

		UnaryOperationAssignment<ComplexDoubleType, DoubleType> imaginaryAssignment = new UnaryOperationAssignment(
				new ComplexImaginaryToRealAdapter());

		int filterIndex = 0;
		for (int a = 0; a < m_halfNumAng; a++) {
			float ang = -a * (float) Math.PI * (1.0f / m_halfNumAng);
			for (int s = 0; s < m_scales.length; s++) {
				for (int f = 0; f < m_frequencies.length; f++) {
					for (int e = 0; e < m_elongations.length; e++) {
						Gabor g = new Gabor(m_radius, ang, m_scales[s], m_frequencies[f], m_elongations[e]);
						for (int evenodd = 0; evenodd < 2; evenodd++) {

							try {
								m_filters[filterIndex] = g.factory().imgFactory(new DoubleType()).create(g,
										new DoubleType());

								if (evenodd == 0) {
									realAssignment.compute(g, m_filters[filterIndex]);
								} else {
									imaginaryAssignment.compute(g, m_filters[filterIndex]);
								}

							} catch (IncompatibleTypeException e1) {
								// TODO handle
								// exception
								e1.printStackTrace();
							}

							if (evenodd == 0) {
								double sub = new Sum<DoubleType, DoubleType>()
										.compute(m_filters[filterIndex].cursor(), new DoubleType()).get()
										/ (g.dimension(0) * g.dimension(1));
								new UnaryConstantRightAssignment(new RealAdd()).compute(m_filters[filterIndex],
										new DoubleType((float) -sub), m_filters[filterIndex]);

								// m_filters[filterIndex]
								// .uAdd(-m_filters[filterIndex].aSum()
								// /
								// (m_filters[filterIndex]
								// .size(0) *
								// m_filters[filterIndex]
								// .size(1)));
							}

							// System.out.println(filterIndex
							// + " "
							// + getFilterIndex(s,
							// f, e, evenodd == 0,
							// a));

							m_names[filterIndex++ % (m_filters.length / m_halfNumAng)] = "Gabor[s=" + m_scales[s]
									+ ";f=" + m_frequencies[f] + ";e=" + m_elongations[e] + ";"
									+ (evenodd == 0 ? "even" : "odd") + "]";

						}
					}

				}
			}

		}

		if (!m_precalcFeatures) {
			m_filterCursors = new Cursor[m_filters.length];
			for (int i = 0; i < m_filters.length; i++) {
				m_filterCursors[i] = m_filters[i].localizingCursor();
			}
		}

		m_numFeatures = m_filters.length / m_halfNumAng;

	}

	// private int getFilterIndex(final int scale, final int freq, final int
	// elon,
	// final boolean even, final int angle) {
	// return angle * m_scales.length * m_frequencies.length
	// * m_elongations.length * 2 + scale * m_frequencies.length
	// * m_elongations.length * 2 + freq * m_elongations.length * 2
	// + elon * 2 + (even ? 0 : 1);
	// }

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void updateImg(Img<T> img) {
		m_img = img;
		m_extendedImg = Views.extendBorder(m_img);
		if (m_precalcFeatures) {
			fftConvolveImage();
		} else {
			m_tmpSrcRandAccess = m_extendedImg.randomAccess();
		}

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void updateDart(Dart dart) {
		m_pos = dart.offset();
		m_currentAngID = dart.directionIndex();

	}

	private final long[] m_kernelRadii = new long[] { m_radius, m_radius };

	/**
	 * {@inheritDoc}
	 */
	@Override
	public double value(int id) {
		double res = -Double.MAX_VALUE; // if (x < m_filterSizeX / 2 ||
										// y
										// <
		// m_filterSizeY / 2
		// || x >= getFeatureTarget().getSizeX() - m_filterSizeX / 2
		// || y >= getFeatureTarget().getSizeY() - m_filterSizeY / 2) {
		// // cursor.getType().set((float)
		// cursor.getType().getMinValue());
		// return res;
		// }
		if (m_precalcFeatures) {
			// res = (float) m_convolvedImgs[m_currentAngID %
			// m_numAng *
			// m_numFeatures
			// + getCurrentFeatureID()].get(x - m_filterSizeX / 2, y
			// - m_filterSizeY / 2);
			RandomAccess<F> tmpFeatRandAccess = m_convolvedImgRandAccess[m_currentAngID % m_halfNumAng * m_numFeatures
					+ id];
			tmpFeatRandAccess.setPosition(m_pos[0], 0);
			tmpFeatRandAccess.setPosition(m_pos[1], 1);
			res = tmpFeatRandAccess.get().getRealDouble();

			// if the current filter is odd and the choosen
			// angle bigger
			// than
			// numAng (only covering a half circle due to
			// the symmetry of
			// the
			// filters)
			if (id % 2 == 1 && m_currentAngID >= m_halfNumAng) {
				res *= -1;
			}

		} else {
			res = convolve(m_tmpSrcRandAccess, m_filterCursors[m_currentAngID % m_halfNumAng * m_numFeatures + id],
					m_pos, m_kernelRadii);

			// if the current filter is odd and the choosen
			// angle bigger
			// than
			// numAng (only covering a half circle due to
			// the symmetry of
			// the
			// filters)
			if (id % 2 == 1 && m_currentAngID >= m_halfNumAng) {
				res *= -1;
			}
		}

		return res;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int numFeatures() {
		return m_names.length;
	}

	// TODO: Is this right?!?!

	@SuppressWarnings("unchecked")
	private void fftConvolveImage() {

		m_convolvedImgs = new Img[m_filters.length];
		m_convolvedImgRandAccess = new RandomAccess[m_filters.length];

		Convolver<T, DoubleType, F> ics = new ImgLib2FourierConvolver<T, DoubleType, F>(m_executorService);

		for (int i = 0; i < m_filters.length; i++) {
			try {
				m_convolvedImgs[i] = m_img.factory().imgFactory(m_featureType).create(m_img, m_featureType);
			} catch (IncompatibleTypeException e) {
			}
			ics.compute(m_extendedImg, m_filters[i], m_convolvedImgs[i]);
			m_convolvedImgRandAccess[i] = m_convolvedImgs[i].randomAccess();
		}
	}

	public void showFilterBank() {

		for (int i = 0; i < m_filters.length; i++) {
			Img<DoubleType> img = m_filters[i];
			int angIdx = i / (m_filters.length / m_numAng);
			double ang = angIdx * Math.PI * (1.0 / m_halfNumAng);
			String name = m_names[i % (m_filters.length / m_halfNumAng)] + "angIdx:" + angIdx + ";ang:" + ang;

			Img<FloatType> floatImg = (Img<FloatType>) KNIPGateway.ops().create().img(img, new FloatType());
			new ImgConvert<DoubleType, FloatType>(new DoubleType(), new FloatType(),
					ImgConvert.ImgConversionTypes.DIRECT, floatImg.factory()).compute(img, floatImg);

			// normalize
			Pair<FloatType, FloatType> minmax = Operations.compute(new MinMax<FloatType>(), floatImg);
			Operations
					.<FloatType, FloatType> map(new Normalize<FloatType>(minmax.getA().getRealDouble(),
							minmax.getB().getRealDouble(), -Float.MAX_VALUE, Float.MAX_VALUE))
					.compute(floatImg, floatImg);

			org.knime.knip.core.awt.AWTImageTools.showInFrame(floatImg, name, 5);

		}

	}

	// /**
	// * {@inheritDoc}
	// */
	// @Override
	// public void readExternal(ObjectInput in) throws IOException,
	// ClassNotFoundException {
	// super.readExternal(in);
	//
	// m_precalcFeatures = in.readBoolean();
	// m_radius = in.readInt();
	// m_numAng = in.readInt();
	//
	// m_scales = new double[in.readInt()];
	// for (int i = 0; i < m_scales.length; i++) {
	// m_scales[i] = in.readDouble();
	// }
	//
	// m_frequencies = new double[in.readInt()];
	// for (int i = 0; i < m_frequencies.length; i++) {
	// m_frequencies[i] = in.readDouble();
	// }
	//
	// m_elongations = new double[in.readInt()];
	// for (int i = 0; i < m_elongations.length; i++) {
	// m_elongations[i] = in.readDouble();
	// }
	// init();
	// }

	// /**
	// * {@inheritDoc}
	// */
	// @Override
	// public void writeExternal(ObjectOutput out) throws IOException {
	// super.writeExternal(out);
	// out.writeBoolean(m_precalcFeatures);
	// out.writeInt(m_radius);
	// out.writeInt(m_numAng);
	//
	// out.writeInt(m_scales.length);
	// for (int i = 0; i < m_scales.length; i++) {
	// out.writeDouble(m_scales[i]);
	// }
	//
	// out.writeInt(m_frequencies.length);
	// for (int i = 0; i < m_frequencies.length; i++) {
	// out.writeDouble(m_frequencies[i]);
	// }
	//
	// out.writeInt(m_elongations.length);
	// for (int i = 0; i < m_elongations.length; i++) {
	// out.writeDouble(m_elongations[i]);
	// }
	//
	// }

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String featureSetId() {
		return "Gabor Filter Feature Factory";
	};

	public static <T extends RealType<T>> void main(final String[] args) {

		// final double[] scales, final double[] frequencies,
		// final double[] elongations, final int radius, final int
		// numAng,
		// boolean precalcFeatures

		// Image<T> img = ImageTools
		// .loadImage("/home/hornm/cell_images/ChHauk-GFP-FAK.tif");
		// Image<DoubleType> fimg = new Convert<T, DoubleType>(new
		// DoubleType(),
		// true)
		// .proc(img);

		GaborFilterFeatureSet<DoubleType, FloatType> gfb = new GaborFilterFeatureSet<DoubleType, FloatType>(
				new double[] { .25f, .5f }, new double[] { 1 }, new double[] { 1 }, 20, 12, true, new FloatType());
		// gfb.getFeatureTarget().setImage(fimg);

		gfb.showFilterBank();
		// gfb.setImage(img);

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void generateFeatImg(RandomAccessibleInterval<F> featImg, int featID, Dart d) {

		int angleID = d.directionIndex();

		Convolver<T, DoubleType, F> ics = new ImgLib2FourierConvolver<T, DoubleType, F>(m_executorService);
		Img<DoubleType> filter = m_filters[angleID % m_halfNumAng * m_numFeatures + featID];

		// if the current filter is odd and the choosen
		// angle bigger
		// than
		// numAng (only covering a half circle due to
		// the symmetry of
		// the
		// filters)
		if (featID % 2 == 1 && angleID >= m_halfNumAng) {
			filter = filter.copy();
			for (DoubleType t : filter) {
				t.set(t.get() * -1);
			}
		}

		ics.compute(m_extendedImg, filter, featImg);

	}

	/**
	 * Straightforward convolution. For small kernels faster than the
	 * convolution in the frequency domain for small images.
	 * 
	 * @param img
	 *            the image in the spatial domain
	 * 
	 * @param kerC
	 *            the kernel in the spatial domain
	 * 
	 * @param pos
	 *            the position to apply the kernel
	 * 
	 * @return
	 */
	public final static <KT extends RealType<KT>, T extends RealType<T>> float convolve(final RandomAccess<T> srcRA,
			final Cursor<KT> kerC, final long[] pos, final long[] kernelRadii) {

		float val = 0;

		srcRA.setPosition(pos);
		kerC.reset();
		while (kerC.hasNext()) {
			kerC.fwd();

			for (int i = 0; i < kernelRadii.length; i++) {
				if (kernelRadii[i] > 0) { // dimension can have
											// zero extension e.g.
											// vertical 1d kernel
					srcRA.setPosition(pos[i] + kerC.getLongPosition(i) - kernelRadii[i], i);
				}
			}

			val += srcRA.get().getRealFloat() * kerC.get().getRealFloat();
		}
		return val;
	}
}
