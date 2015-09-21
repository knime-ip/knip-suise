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
package org.knime.knip.suise.node.boundarymodel;

import java.awt.Image;
import java.awt.Polygon;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.DefaultLinearAxis;
import net.imagej.space.DefaultCalibratedSpace;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgView;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.logic.BoolType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.util.ConstantUtils;

import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.RowIterator;
import org.knime.core.data.collection.CollectionCellFactory;
import org.knime.core.data.collection.ListCell;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.PortObjectSpec;
import org.knime.core.node.port.PortType;
import org.knime.knip.base.data.img.ImgPlusCell;
import org.knime.knip.base.data.img.ImgPlusCellFactory;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.data.labeling.LabelingValue;
import org.knime.knip.base.exceptions.KNIPException;
import org.knime.knip.base.node.NodeUtils;
import org.knime.knip.base.node.TwoValuesToCellNodeModel;
import org.knime.knip.core.KNIPGateway;
import org.knime.knip.core.data.img.DefaultImgMetadata;
import org.knime.knip.core.util.PolygonTools;
/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import org.knime.knip.suise.node.pixclassmodel.SettingsModelWekaClassifier;
import org.knime.knip.suise.node.pixfeat2d.PixFeaturesNodeModel;
import org.knime.knip.suise.node.pixfeat2d.angledep.BufferedPixelFeatureSet;

import weka.classifiers.trees.RandomForest;

/**
 * 
 * @author hornm, University of Konstanz
 */
public class BoundaryModelNodeModel<F extends RealType<F>, T extends RealType<T>, L extends Comparable<L>>
		extends TwoValuesToCellNodeModel<ImgPlusValue<F>, LabelingValue<L>, ListCell> {

	private static final NodeLogger LOGGER = NodeLogger.getLogger(BoundaryModelNodeModel.class);

	private static final int DIM_X = 0;

	private static final int DIM_Y = 1;

	static SettingsModelString createSrcImgColModel() {
		return new SettingsModelString("src_img_column", "");
	}

	static SettingsModelString createApplyImgColModel() {
		return new SettingsModelString("apply_img_column", "");
	}

	static SettingsModelString createMaskImgColModel() {
		return new SettingsModelString("mask_img_column", "");
	}

	static SettingsModelString createSelectionStrategyModel() {
		return new SettingsModelString("selection_strategy",
				BoundaryModel.SelectionStrategy.ITERATIVE_INFERENCE_SELECTION.name());
	}

	static SettingsModelWekaClassifier createClassifierModel() {
		return new SettingsModelWekaClassifier("weka_classifier", new RandomForest());
	}

	static SettingsModelString createOptionalParametersModel() {
		return new SettingsModelString("optional_parameters",
				BoundaryModel.OPTIONAL_PARAMETER_STDEV + "=200," + BoundaryModel.OPTIONAL_PARAMETER_BIAS + "=100");
	}

	final static SettingsModelString createOrientationDimLabelModel() {
		return new SettingsModelString("orientation_dim_label", PixFeaturesNodeModel.DEFAULT_ORIENTATION_DIM_LABEL);
	}

	final static SettingsModelString createFeatureDimLabelModel() {
		return new SettingsModelString("feature_dim_label", PixFeaturesNodeModel.DEFAULT_FEATURE_DIM_LABEL);
	}

	private SettingsModelString m_srcImgCol = createSrcImgColModel();

	private SettingsModelString m_applyImgCol = createApplyImgColModel();

	private SettingsModelString m_maskImgCol = createMaskImgColModel();

	private SettingsModelString m_selectionStrategy = createSelectionStrategyModel();

	private SettingsModelWekaClassifier m_classifier = createClassifierModel();

	private SettingsModelString m_optionalParameters = createOptionalParametersModel();

	private SettingsModelString m_orientationDimLabel = createOrientationDimLabelModel();

	private SettingsModelString m_featureDimLabel = createFeatureDimLabelModel();

	private BoundaryModel<T, F> m_boundaryModel;

	private ImgPlusCellFactory m_imgCellFactory;

	public BoundaryModelNodeModel() {
		super(new PortType[] { new PortType(BufferedDataTable.class, true) });
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void addSettingsModels(List<SettingsModel> settingsModels) {
		settingsModels.add(m_srcImgCol);
		settingsModels.add(m_applyImgCol);
		settingsModels.add(m_maskImgCol);
		settingsModels.add(m_selectionStrategy);
		settingsModels.add(m_classifier);
		settingsModels.add(m_optionalParameters);
		settingsModels.add(m_orientationDimLabel);
		settingsModels.add(m_featureDimLabel);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected PortObjectSpec[] configure(PortObjectSpec[] inSpecs) throws InvalidSettingsException {
		if (inSpecs[1] == null) {
			return super.configure(inSpecs);
		} else {
			getApplyImgColIndex((DataTableSpec) inSpecs[1]);
			return new PortObjectSpec[] { new DataTableSpec(
					new DataColumnSpecCreator("Prob. Map", ListCell.getCollectionType(ImgPlusCell.TYPE))
							.createSpec()) };
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected PortObject[] execute(PortObject[] inObjects, ExecutionContext exec) throws Exception {
		exec.setProgress("Build boundary model ...");

		BufferedDataTable inTable0 = (BufferedDataTable) inObjects[0];

		RowIterator rowIterator = inTable0.iterator();

		int firstColIdx = getFirstColumnIdx(inTable0.getDataTableSpec());
		int secondColidx = getSecondColumnIdx(inTable0.getDataTableSpec());

		int srcImgColIndex = getSrcImgColIndex(inTable0.getDataTableSpec());

		m_boundaryModel = new BoundaryModel<T, F>();

		// set selection strategy
		m_boundaryModel
				.setSelectionStrategy(BoundaryModel.SelectionStrategy.valueOf(m_selectionStrategy.getStringValue()));

		// set classifier
		m_boundaryModel.setWekaClassifier(m_classifier.getClassifier());

		// set optional parameters
		String[] splits = m_optionalParameters.getStringValue().split(",");
		HashMap<String, String> map = new HashMap<String, String>(splits.length);
		for (String s : splits) {
			String[] keyvalue = s.split("=");
			if (keyvalue.length == 2) {
				map.put(keyvalue[0].trim(), keyvalue[1].trim());
			}
		}
		m_boundaryModel.setOptionalParameters(map);

		while (rowIterator.hasNext()) {
			DataRow row = rowIterator.next();

			if (row.getCell(firstColIdx).isMissing() || row.getCell(secondColidx).isMissing()) {
				LOGGER.warn("Missing cells in row " + row.getKey() + ". Row has been skipped.");
				continue;
			}

			// the feature image
			@SuppressWarnings("unchecked")
			ImgPlus<F> img = ((ImgPlusValue<F>) row.getCell(firstColIdx)).getImgPlus();
			Img<T> srcImg = null;
			if (srcImgColIndex > -1) {
				srcImg = ((ImgPlusValue<T>) row.getCell(srcImgColIndex)).getImgPlus();
				if (srcImg.numDimensions() != 2) {
					LOGGER.warn("The source image for debugging in row " + row.getKey()
							+ " has been skipped as its not 2-dimensional!");
					srcImg = null;
				}
			}

			int dimFeat = img.dimensionIndex(Axes.get(m_featureDimLabel.getStringValue()));
			int dimAngle = img.dimensionIndex(Axes.get(m_orientationDimLabel.getStringValue()));
			if (dimFeat == -1 || dimAngle == -1) {
				throw new KNIPException(
						"At least one of the dimension label (orientation/feature) is not present in the image "
								+ img.getName() + ".");
			}

			BufferedPixelFeatureSet<F> featSet = new BufferedPixelFeatureSet<F>(
					new String[(int) img.dimension(dimFeat)], DIM_X, DIM_Y, dimFeat, dimAngle);
			m_boundaryModel.setFeatureSet(featSet);

			LabelingValue<L> labVal = (LabelingValue<L>) row.getCell(secondColidx);

			featSet.updateImg(img);

			// generate polygons and add them to the boundary model
			for (Polygon p : getPolygons(labVal.getLabeling())) {
				m_boundaryModel.addSamples(p, (int) img.dimension(0), (int) img.dimension(1), srcImg);
			}

		}
		// building the cell boundary model
		m_boundaryModel.buildModel();

		exec.setMessage("Classifying pixels ...");

		m_imgCellFactory = new ImgPlusCellFactory(exec);

		if (inObjects[1] == null) {
			// if the optional inport is not connected than classify the
			// training examples
			return super.execute(inObjects, exec);
		} else {
			// else classify the so far unseen images from the optional inport
			BufferedDataTable inTable1 = (BufferedDataTable) inObjects[1];
			BufferedDataContainer dataCon = exec.createDataContainer(new DataTableSpec(
					new DataColumnSpecCreator("Prob. Map", ListCell.getCollectionType(ImgPlusCell.TYPE)).createSpec()));
			RowIterator rowIt = inTable1.iterator();
			int colidx = getApplyImgColIndex(inTable1.getDataTableSpec());
			int maskImgColIndex = getMaskImgColIndex(inTable1.getDataTableSpec());
			Img<BitType> mask = null;
			int rowCount = 0;
			while (rowIt.hasNext()) {
				DataRow row = rowIt.next();
				if (maskImgColIndex != -1) {
					mask = ((ImgPlusValue<BitType>) row.getCell(maskImgColIndex)).getImgPlus();
				}
				DefaultRow res = new DefaultRow(row.getKey(),
						computeResult((ImgPlusValue<F>) row.getCell(colidx), mask));
				dataCon.addRowToTable(res);
				exec.checkCanceled();
				exec.setProgress((double) rowCount++ / inTable1.getRowCount());

			}
			dataCon.close();
			return new BufferedDataTable[] { dataCon.getTable() };
		}
	}

	private int getSrcImgColIndex(DataTableSpec inSpec) throws InvalidSettingsException {
		int colIdx = -1;
		if (m_srcImgCol.getStringValue() != null) {
			colIdx = NodeUtils.autoColumnSelection(inSpec, m_srcImgCol, ImgPlusValue.class, this.getClass());
		}
		return colIdx;
	}

	private int getApplyImgColIndex(DataTableSpec inSpec) throws InvalidSettingsException {
		int colIdx = -1;
		if (m_applyImgCol.getStringValue() != null) {
			colIdx = NodeUtils.autoColumnSelection(inSpec, m_applyImgCol, ImgPlusValue.class, this.getClass());
		}
		return colIdx;
	}

	private int getMaskImgColIndex(DataTableSpec inSpec) throws InvalidSettingsException {
		int colIdx = -1;
		if (m_maskImgCol.getStringValue() != null) {
			colIdx = NodeUtils.autoColumnSelection(inSpec, m_maskImgCol, ImgPlusValue.class, this.getClass());
		}
		return colIdx;
	}

	/* Creates a list of polygons from a labeling */
	private Collection<Polygon> getPolygons(RandomAccessibleInterval<LabelingType<L>> randomAccessibleInterval) {
		int[] offset = new int[2];
		long[] tmp = new long[2];
		ArrayList<Polygon> polygons = new ArrayList<Polygon>();
		final LabelRegions<L> regions = KNIPGateway.regions().regions(randomAccessibleInterval);
		for (final L label : regions.getExistingLabels()) {
			LabelRegion<L> region = regions.getLabelRegion(label);
			Img<BitType> mask = binaryMask(region);
			region.min(tmp);
			offset[0] = (int) tmp[0];
			offset[1] = (int) tmp[1];
			Polygon p = PolygonTools.extractPolygon(mask, offset);
			polygons.add(p);
		}
		return polygons;
	}

	/* Creates the binary mask from an iterable interval */
	private Img<BitType> binaryMask(IterableInterval<Void> ii) {
		Img<BitType> binaryMask = new ArrayImgFactory<BitType>().create(ii, new BitType());
		RandomAccess<BitType> maskRA = binaryMask.randomAccess();

		Cursor<Void> cur = ii.localizingCursor();
		while (cur.hasNext()) {
			cur.fwd();
			for (int d = 0; d < cur.numDimensions(); d++) {
				maskRA.setPosition(cur.getLongPosition(d) - ii.min(d), d);
			}
			maskRA.get().set(true);

		}
		return binaryMask;

	}

	protected Image getBoundaryModelCreationDebugImage() {
		if (m_boundaryModel != null) {
			return m_boundaryModel.getDebugImage();
		} else {
			return null;
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected ListCell compute(ImgPlusValue<F> cellValue1, LabelingValue<L> cellValue2) throws Exception {
		/*
		 * Apply the boundary classifier to the source images and/or the cell
		 * samples, if desired
		 */
		return computeResult(cellValue1, null);

	}

	private ListCell computeResult(ImgPlusValue<F> val, Img<BitType> mask) throws Exception {
		// the source feature image
		ImgPlus<F> img = val.getImgPlus();
		int dimFeat = img.dimensionIndex(Axes.get(m_featureDimLabel.getStringValue()));
		int dimAngle = img.dimensionIndex(Axes.get(m_orientationDimLabel.getStringValue()));
		if (dimFeat == -1 || dimAngle == -1) {
			throw new KNIPException("One of the dimension labels (orientation/feature) are not present in the image "
					+ img.getName() + ".");
		}
		BufferedPixelFeatureSet<F> featSet = new BufferedPixelFeatureSet<F>(new String[(int) img.dimension(dimFeat)],
				DIM_X, DIM_Y, dimFeat, dimAngle);
		featSet.updateImg(img);
		m_boundaryModel.setFeatureSet(featSet);

		if (mask == null) {
			mask =

			new ImgView<BitType>(ConstantUtils.constantRandomAccessibleInterval(new BitType(true), 2,
					new FinalInterval(new long[] { img.dimension(0), img.dimension(1) })), null);
		}
		Img<ByteType>[] classResultImg = m_boundaryModel.classifyImageContourModelwise(new int[2], mask);

		// create cells
		ArrayList<ImgPlusCell> cells = new ArrayList<ImgPlusCell>(classResultImg.length);
		for (int i = 0; i < classResultImg.length; i++) {
			ImgPlusCell<ByteType> cell = m_imgCellFactory
					.createCell(
							new ImgPlus<>(classResultImg[i],
									new DefaultImgMetadata(
											new DefaultCalibratedSpace(new DefaultLinearAxis(Axes.X),
													new DefaultLinearAxis(Axes.Y),
													new DefaultLinearAxis(
															Axes.get(m_orientationDimLabel.getStringValue()))),
											img, img, img)));
			cells.add(cell);
		}

		return CollectionCellFactory.createListCell(cells);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataType getOutDataCellListCellType() {
		return ImgPlusCell.TYPE;
	}

}
