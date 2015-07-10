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
package org.knime.knip.suise.node.pixfeat2d;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.imagej.ImgPlus;
import net.imagej.ImgPlusMetadata;
import net.imagej.axis.Axes;
import net.imagej.axis.CalibratedAxis;
import net.imagej.axis.DefaultLinearAxis;
import net.imagej.space.DefaultCalibratedSpace;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.ops.operation.SubsetOperations;
import net.imglib2.ops.operation.img.unary.ImgConvert;
import net.imglib2.ops.operation.img.unary.ImgConvert.ImgConversionTypes;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;

import org.knime.core.node.ExecutionContext;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;
import org.knime.knip.base.data.img.ImgPlusCell;
import org.knime.knip.base.data.img.ImgPlusCellFactory;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.node.ValueToCellNodeModel;
import org.knime.knip.core.KNIPGateway;
import org.knime.knip.core.data.img.DefaultImgMetadata;
import org.knime.knip.core.ops.img.ImgPlusNormalize;
import org.knime.knip.suise.data.feat.PixFeatureFactory;
import org.knime.knip.suise.data.feat.PixFeatureSet;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class PixFeaturesNodeModel<T extends RealType<T>> extends
		ValueToCellNodeModel<ImgPlusValue<T>, ImgPlusCell<UnsignedByteType>> {

	public static final String DEFAULT_FEATURE_DIM_LABEL = "F";

	public static final String DEFAULT_ORIENTATION_DIM_LABEL = "O";

	final static SettingsModelStringArray createActiveFeatureSetModel() {
		return new SettingsModelStringArray("active featureSets", new String[0]);
	}

	final static SettingsModelIntegerBounded createNumOrientationsModel(
			int maxNumOrientations) {
		return new SettingsModelIntegerBounded("max_num_orientations",
				maxNumOrientations, 1, maxNumOrientations);
	}

	final static SettingsModelString createOrientationDimLabelModel() {
		return new SettingsModelString("orientation_dim_label",
				DEFAULT_ORIENTATION_DIM_LABEL);
	}

	final static SettingsModelString createFeatureDimLabelModel() {
		return new SettingsModelString("feature_dim_label",
				DEFAULT_FEATURE_DIM_LABEL);
	}

	private SettingsModelStringArray m_activeFeatureSets = createActiveFeatureSetModel();

	private SettingsModelIntegerBounded m_smNumOrientations;

	private SettingsModelString m_smFeatDimLabel;

	private SettingsModelString m_smOrientationDimLabel;

	private PixFeatureFactory<T> m_featFac;

	private final PixFeatureSetProvider[] m_pixFeatProviders;

	private final int m_pixFeatDimensionality;

	private ImgPlusCellFactory m_imgCellFactory;

	public PixFeaturesNodeModel(PixFeatureSetProvider[] pixFeatProviders,
			int pixFeatDim, int maxNumOrientations) {
		super();
		m_pixFeatProviders = pixFeatProviders;
		m_pixFeatDimensionality = pixFeatDim;
		m_smNumOrientations = createNumOrientationsModel(maxNumOrientations);
		m_smFeatDimLabel = createFeatureDimLabelModel();
		m_smOrientationDimLabel = createOrientationDimLabelModel();
		if (maxNumOrientations == 1) {
			m_smNumOrientations.setEnabled(false);
			m_smOrientationDimLabel.setEnabled(false);
		}

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void addSettingsModels(List<SettingsModel> settingsModels) {
		settingsModels.add(m_activeFeatureSets);
		settingsModels.add(m_smNumOrientations);
		settingsModels.add(m_smOrientationDimLabel);
		settingsModels.add(m_smFeatDimLabel);
		for (PixFeatureSetProvider<T> p : m_pixFeatProviders) {
			p.initAndAddSettingsModels(settingsModels);
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void prepareExecute(ExecutionContext exec) {
		// create feature factory
		Set<String> active = new HashSet<String>(
				Arrays.asList(m_activeFeatureSets.getStringArrayValue()));
		List<PixFeatureSet<T>> featSets = new ArrayList<PixFeatureSet<T>>(
				m_pixFeatProviders.length);
		for (PixFeatureSetProvider<T> p : m_pixFeatProviders) {
			if (active.contains(p.getFeatureSetId())) {
				featSets.add(p.getPixFeatureSet(m_smNumOrientations
						.getIntValue()));
			}
		}

		m_featFac = new PixFeatureFactory<T>(featSets);

		m_imgCellFactory = new ImgPlusCellFactory(exec);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected ImgPlusCell<UnsignedByteType> compute(ImgPlusValue<T> cellValue)
			throws Exception {

		ImgPlus<FloatType> featImg = new ImgPlus(m_featFac.makeFeatureImage(
				cellValue.getImgPlus(), m_pixFeatDimensionality,
				m_smNumOrientations.getIntValue()));

		// normalize image interval plane-wise
		SubsetOperations.iterate(new ImgPlusNormalize<FloatType>(0,
				new FloatType(), null, false), new int[] { 0, 1 }, featImg,
				featImg);

		Img<UnsignedByteType> res = (Img<UnsignedByteType>) KNIPGateway.ops()
				.createImg(featImg, new UnsignedByteType());
		ImgConvert<FloatType, UnsignedByteType> convert = new ImgConvert<FloatType, UnsignedByteType>(
				new FloatType(), new UnsignedByteType(),
				ImgConversionTypes.SCALE, res.factory());
		convert.compute((RandomAccessibleInterval<FloatType>) featImg, res);
		ImgPlusMetadata metadata = cellValue.getMetadata();
		if (m_smNumOrientations.getIntValue() > 1) {
			List<CalibratedAxis> axes = new ArrayList<CalibratedAxis>(
					metadata.numDimensions() + 1);
			for (int i = 0; i < metadata.numDimensions(); i++) {
				axes.add(metadata.axis(i));
			}
			axes.add(new DefaultLinearAxis(Axes.get(m_smFeatDimLabel
					.getStringValue())));
			axes.add(new DefaultLinearAxis(Axes.get(m_smOrientationDimLabel
					.getStringValue())));
			return m_imgCellFactory.createCell(res, new DefaultImgMetadata(
					new DefaultCalibratedSpace(axes), metadata, metadata,
					metadata));
		} else {
			List<CalibratedAxis> axes = new ArrayList<CalibratedAxis>(
					metadata.numDimensions() + 1);
			for (int i = 0; i < metadata.numDimensions(); i++) {
				axes.add(metadata.axis(i));
			}
			axes.add(new DefaultLinearAxis(Axes.get(m_smFeatDimLabel
					.getStringValue())));
			return m_imgCellFactory.createCell(res, new DefaultImgMetadata(
					new DefaultCalibratedSpace(axes), metadata, metadata,
					metadata));
		}
	}
}
