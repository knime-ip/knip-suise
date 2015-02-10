/*
 * ------------------------------------------import java.util.List;

import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.converter.Converter;
import net.imglib2.converter.Converters;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.histogram.Integer1dBinMapper;
import net.imglib2.img.Img;
import net.imglib2.img.ImgView;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.meta.ImgPlus;
import net.imglib2.ops.operation.iterableinterval.unary.MakeHistogram;
import net.imglib2.ops.operation.randomaccessible.binary.FloodFill;
import net.imglib2.ops.operation.randomaccessible.unary.FillHoles;
import net.imglib2.ops.operation.randomaccessibleinterval.unary.regiongrowing.CCA;
import net.imglib2.ops.types.ConnectedType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

import org.knime.core.data.DataCell;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.data.labeling.LabelingCell;
import org.knime.knip.base.data.labeling.LabelingCellFactory;
import org.knime.knip.base.exceptions.ImageTypeNotCompatibleException;
import org.knime.knip.base.exceptions.KNIPRuntimeException;
import org.knime.knip.base.node.ValueToCellNodeDialog;
import org.knime.knip.base.node.ValueToCellNodeFactory;
import org.knime.knip.base.node.ValueToCellNodeModel;
import org.knime.knip.core.data.img.DefaultLabelingMetadata;
import org.knime.knip.core.util.NeighborhoodUtils;
 are deemed to be separate and independent programs and to not be
 *  covered works.  Notwithstanding anything to the contrary in the
 *  License, the License does not apply to Nodes, you are not required to
 *  license Nodes under the License, and you are granted a license to
 *  prepare and propagate Nodes, in each case even if such Nodes are
 *  propagated with or for interoperation with KNIME. The owner of a Node
 *  may freely choose the license terms applicable to such Node, including
 *  when such Node is propagated with or for interoperation with KNIME.
 * ------------------------------------------------------------------------
 * 
 * History
 *   Jun 3, 2013 (hornm): created
 */

package org.knime.knip.suise.node.ucmregions;

import java.util.List;

import net.imagej.ImgPlus;
import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.converter.Converter;
import net.imglib2.converter.Converters;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.histogram.Integer1dBinMapper;
import net.imglib2.img.Img;
import net.imglib2.img.ImgView;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.ops.operation.iterableinterval.unary.MakeHistogram;
import net.imglib2.ops.operation.randomaccessible.binary.FloodFill;
import net.imglib2.ops.operation.randomaccessible.unary.FillHoles;
import net.imglib2.ops.operation.randomaccessibleinterval.unary.regiongrowing.CCA;
import net.imglib2.ops.types.ConnectedType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

import org.knime.core.data.DataCell;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.data.labeling.LabelingCell;
import org.knime.knip.base.data.labeling.LabelingCellFactory;
import org.knime.knip.base.exceptions.ImageTypeNotCompatibleException;
import org.knime.knip.base.exceptions.KNIPRuntimeException;
import org.knime.knip.base.node.ValueToCellNodeDialog;
import org.knime.knip.base.node.ValueToCellNodeFactory;
import org.knime.knip.base.node.ValueToCellNodeModel;
import org.knime.knip.core.data.img.DefaultLabelingMetadata;
import org.knime.knip.core.util.NeighborhoodUtils;

/**
 * Extracts all possible regions from an ultra-metric contour map (UCM) by
 * thresholding at each available intensity level.
 * 
 * @author hornm, University of Konstanz
 */
public class UCMRegionsNodeFactory<T extends IntegerType<T>> extends
		ValueToCellNodeFactory<ImgPlusValue<T>> {

	private static final SettingsModelInteger createMinBoundaryPixNumberModel() {
		return new SettingsModelInteger("minimum_boundary_pix_number", 0);
	}

	private static final SettingsModelInteger createMinPixIntensityModel() {
		return new SettingsModelInteger("minimum_pixel_intensity", 0);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected ValueToCellNodeDialog<ImgPlusValue<T>> createNodeDialog() {
		return new ValueToCellNodeDialog<ImgPlusValue<T>>() {

			@Override
			public void addDialogComponents() {
				addDialogComponent("Options", "", new DialogComponentNumber(
						createMinBoundaryPixNumberModel(),
						"minimum number of boundary pixels", 1));
				addDialogComponent("Options", "", new DialogComponentNumber(
						createMinPixIntensityModel(),
						"minimum pixel intensity", 1));

			}
		};
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public ValueToCellNodeModel<ImgPlusValue<T>, ? extends DataCell> createNodeModel() {
		return new ValueToCellNodeModel<ImgPlusValue<T>, LabelingCell<Integer>>() {

			private SettingsModelInteger m_smMinBoundaryPixNum = createMinBoundaryPixNumberModel();

			private SettingsModelInteger m_smMinPixIntensity = createMinPixIntensityModel();

			private LabelingCellFactory m_cellFactory;

			private ExecutionContext m_exec;

			@Override
			protected void addSettingsModels(final List<SettingsModel> settingsModels) {
				settingsModels.add(m_smMinBoundaryPixNum);
				settingsModels.add(m_smMinPixIntensity);

			}

			/**
			 * {@inheritDoc}
			 */
			@Override
			protected void prepareExecute(final ExecutionContext exec) {
				m_cellFactory = new LabelingCellFactory(exec);
				m_exec = exec;
			}

			@Override
			protected LabelingCell<Integer> compute(final ImgPlusValue<T> cellValue)
					throws Exception {
				final ImgPlus<T> img = cellValue.getImgPlus();
				if (!(img.firstElement().createVariable() instanceof IntegerType)) {
					throw new ImageTypeNotCompatibleException("UCMRegions",
							img.firstElement(), IntegerType.class);
				}
				final Img<BitType> levelComps = img.factory()
						.imgFactory(new BitType()).create(img, new BitType());
				final Img<BitType> tmp = levelComps.factory().create(img,
						new BitType());

				Img<BitType> th0 = null;
				Img<BitType> th1 = null;

				final Histogram1d<T> hist = new Histogram1d<T>(
						new Integer1dBinMapper<T>((int)img.firstElement().getMinValue(), (int)img.firstElement().getMaxValue(), false));
				new MakeHistogram<T>((int) hist.getBinCount()).compute(img,
						hist);

				final FloodFill<BitType> floodFill = new FloodFill<BitType>(
						ConnectedType.EIGHT_CONNECTED);
				final FillHoles fillHoles = new FillHoles(
						ConnectedType.FOUR_CONNECTED);
				final CCA<BitType> cca = new CCA<BitType>(
						NeighborhoodUtils.get4ConStructuringElement(img
								.numDimensions()), new BitType());
				final Labeling<Integer> res = new NativeImgLabeling<Integer, IntType>(
						img.factory().imgFactory(new IntType())
								.create(img, new IntType()));
				for (long i = hist.getBinCount() - 1; i >= m_smMinPixIntensity
						.getIntValue(); i--) {
					if (hist.frequency(i) <= m_smMinBoundaryPixNum
							.getIntValue()) {
						continue;
					}
					
					m_exec.checkCanceled();
//					break;

					System.out.println(hist.frequency(i));
					
					final ThresholdConverter<T> c = new ThresholdConverter<T>(i);

					// thresholded image at all available levels
					th1 = new ImgView<BitType>(
							Converters.convert(
									(RandomAccessibleInterval<T>) img, c,
									new BitType()), img.factory().imgFactory(
									new BitType()));
					final Cursor<BitType> th1C = th1.localizingCursor();
					final Cursor<BitType> th0C = th0 == null ? null : th0.cursor();
					Cursor<BitType> lcC = levelComps.cursor();
					while (th1C.hasNext()) {
						th1C.fwd();
						lcC.fwd();
						if (th0C != null) {
							th0C.fwd();
						}
						if (th1C.get().get()
								&& (th0C == null || !th0C.get().get())
								&& !lcC.get().get()) {
							// if the pixel was newly added at the
							// current thresholding level -> use position for
							// flood filling on the current level
							floodFill.compute(th1, th1C, levelComps);

						}

					}
					th0 = th1;
					// fill holes on the components at the current level
					fillHoles.compute(levelComps, tmp);
					// remove the boundaries
					final Cursor<BitType> tmpC = tmp.cursor();
					lcC = levelComps.cursor();
					while (lcC.hasNext()) {
						lcC.fwd();
						tmpC.fwd();
						if (tmpC.get().get() && lcC.get().get()) {
							tmpC.get().set(false);
						}
					}

					// perform the connected component analysis
					cca.compute(tmp, res);

					for (final BitType t : levelComps) {
						t.set(false);
					}
				}
				if (res.getLabels().size() == 0) {
					throw new KNIPRuntimeException("Empty labeling!");
				} else {
					return m_cellFactory.createCell(res,
							new DefaultLabelingMetadata(img, img, img, null));
				}
			}
		};
	}

	private class ThresholdConverter<T extends RealType<T>> implements
			Converter<T, BitType> {

		private double m_threshold;

		private BitType m_type;

		public ThresholdConverter(final double threshold) {
			m_threshold = threshold;
			m_type = new BitType();

		}

		public void setThreshold(final double t) {
			m_threshold = t;
		}

		/**
		 * {@inheritDoc}
		 */
		@Override
		public void convert(final T input, final BitType output) {
			output.set(input.getRealDouble() >= m_threshold);

		}
	}

}
