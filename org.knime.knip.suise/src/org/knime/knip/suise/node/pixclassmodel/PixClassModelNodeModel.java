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
package org.knime.knip.suise.node.pixclassmodel;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.imagej.ImgPlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.labeling.Labeling;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Util;

import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowIterator;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.PortObjectSpec;
import org.knime.core.node.port.PortType;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.data.labeling.LabelingValue;
import org.knime.knip.base.exceptions.KNIPRuntimeException;
import org.knime.knip.base.node.NodeUtils;
import org.knime.knip.base.node.nodesettings.SettingsModelDimSelection;
import org.knime.knip.core.KNIPGateway;
import org.knime.knip.suise.node.port.WekaClassifierPortObject;
import org.knime.knip.suise.node.port.WekaClassifierPortObjectSpec;
import org.knime.knip.suise.ops.BuildTrainingData;

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.Classifier;
import weka.classifiers.trees.RandomForest;
import weka.core.Instance;
import weka.core.Instances;

/**
 * 
 * @author Martin Horn, University of Konstanz
 */
public class PixClassModelNodeModel<T extends RealType<T>, L extends Comparable<L>>
		extends NodeModel {

	protected static final NodeLogger LOGGER = NodeLogger
			.getLogger(PixClassModelNodeModel.class);

	static final SettingsModelString createLabColumnModel() {
		return new SettingsModelString("labeling_column", "");
	}

	static final SettingsModelString createImgColumnModel() {
		return new SettingsModelString("img_column", "");
	}

	static final SettingsModelDimSelection createDimSelectionModel() {
		return new SettingsModelDimSelection("dim_selection", "X", "Y");
	}

	static final SettingsModelDimSelection createFeatDimSelectionModel() {
		return new SettingsModelDimSelection("feat_dim_selection", "C");
	}

	static final SettingsModelWekaClassifier createClassifierSelectionModel() {
		return new SettingsModelWekaClassifier("weka_classifier",
				new RandomForest());
	}

	static final SettingsModelDouble createSampleRateModel() {
		return new SettingsModelDouble("sample_rate", 0.1);
	}

	static final SettingsModelBoolean createBalanceModel() {
		return new SettingsModelBoolean("balance_instances", false);
	}

	public static final int INPUT_DATATABLE_PORT = 0;

	private SettingsModelString m_labelingColumn = createLabColumnModel();

	private SettingsModelString m_imgColumn = createImgColumnModel();

	private SettingsModelDimSelection m_dimSelection = createDimSelectionModel();

	private SettingsModelDimSelection m_featDimSelection = createFeatDimSelectionModel();

	private SettingsModelWekaClassifier m_classifierSelection = createClassifierSelectionModel();

	private SettingsModelDouble m_resampleRate = createSampleRateModel();

	private SettingsModelBoolean m_balanceClassInstances = createBalanceModel();

	protected PixClassModelNodeModel() {
		super(new PortType[] { BufferedDataTable.TYPE },
				new PortType[] { WekaClassifierPortObject.TYPE });
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected PortObjectSpec[] configure(PortObjectSpec[] inSpecs)
			throws InvalidSettingsException {

		getImgColumnIndex((DataTableSpec) inSpecs[0]);
		getLabelingColumnIndex((DataTableSpec) inSpecs[0]);
		return null;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected PortObject[] execute(PortObject[] inObjects, ExecutionContext exec)
			throws Exception {

		BufferedDataTable inTable = (BufferedDataTable) inObjects[0];

		int imgColIdx = getImgColumnIndex(inTable.getDataTableSpec());
		int labColIdx = getLabelingColumnIndex(inTable.getDataTableSpec());

		// retrieve all available labels
		RowIterator it = inTable.iterator();
		DataRow row;
		Set<L> labels = new HashSet<L>();
		Instances trainingSet = null;
		int rowCount = inTable.getRowCount();
		int i = 0;
		while (it.hasNext()) {
			row = it.next();
			if (row.getCell(labColIdx).isMissing()
					|| row.getCell(imgColIdx).isMissing()) {
				setWarningMessage("Errors occurred while execution! See console for details.");
				LOGGER.warn("Missing cell in row " + row.getKey()
						+ ". Row skipped!");
				continue;
			}
			RandomAccessibleInterval<LabelingType<L>> lab = ((LabelingValue<L>) row
					.getCell(labColIdx)).getLabeling();
			ImgPlus<T> img = ((ImgPlusValue<T>) row.getCell(imgColIdx))
					.getImgPlus();

			// collect available labels
			LabelRegions<L> regions = KNIPGateway.regions().regions(lab);
			labels.addAll(regions.getExistingLabels());

			int[] tmp = m_featDimSelection.getSelectedDimIndices(
					img.numDimensions(), img);
			if (tmp.length == 0) {
				setWarningMessage("Errors occurred while execution! See console for details.");
				LOGGER.warn("Feature dimensions doesn't exist in image in row "
						+ row.getKey() + ". Row skipped!");
				continue;
			}
			int featDim = tmp[0];

			int[] dimIndices = m_dimSelection.getSelectedDimIndices(
					img.numDimensions(), img);
			List<String> classLabels = new ArrayList<String>();
			for (L label : regions.getExistingLabels()) {
				classLabels.add(label.toString());
			}
			BuildTrainingData<L, T> btd = new BuildTrainingData<L, T>(
					classLabels, dimIndices, featDim,
					m_resampleRate.getDoubleValue(),
					m_balanceClassInstances.getBooleanValue());

			if (trainingSet == null) {
				trainingSet = btd.bufferFactory().instantiate(lab, img);
			}
			exec.setProgress("Building training set for row " + row.getKey());
			try {
				btd.compute(lab, img, trainingSet);
			} catch (KNIPRuntimeException e) {
				setWarningMessage("Errors occurred while execution! See console for details.");
				LOGGER.warn("Row " + row.getKey() + " skipped. "
						+ e.getLocalizedMessage());
			}

			exec.checkCanceled();
			exec.setProgress((double) i / rowCount);
			i++;
		}

		// build classifier
		exec.setProgress("Build classifier ...");
		if (trainingSet == null) {
			throw new IllegalStateException(
					"No training set could be created due to the lack of training samples. Maybe wrong (i.e. non-existent) feature dimension selected!?");
		}

		// count instances per class for debugging purposes
		double[] classDistr = new double[trainingSet.numClasses()];
		for (Instance instance : trainingSet) {
			classDistr[(int) instance.classValue()]++;
		}
		Classifier classifier = m_classifierSelection.getClassifier();
		classifier.buildClassifier(trainingSet);
		return new PortObject[] { new WekaClassifierPortObject(classifier,
				trainingSet, new WekaClassifierPortObjectSpec(
						labels.toArray(new String[labels.size()]))) };

	}

	/**
	 * @param inSpec
	 * @return index of the image column
	 * @throws InvalidSettingsException
	 */
	protected final int getImgColumnIndex(DataTableSpec inSpec)
			throws InvalidSettingsException {
		return NodeUtils.silentOptionalAutoColumnSelection(inSpec, m_imgColumn,
				ImgPlusValue.class);

	}

	/**
	 * @param inSpec
	 * @return index of the labeling column
	 * @throws InvalidSettingsException
	 */
	protected final int getLabelingColumnIndex(DataTableSpec inSpec)
			throws InvalidSettingsException {
		return NodeUtils.silentOptionalAutoColumnSelection(inSpec,
				m_labelingColumn, LabelingValue.class);

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(NodeSettingsWO settings) {
		m_imgColumn.saveSettingsTo(settings);
		m_labelingColumn.saveSettingsTo(settings);
		m_dimSelection.saveSettingsTo(settings);
		m_featDimSelection.saveSettingsTo(settings);
		m_classifierSelection.saveSettingsTo(settings);
		m_resampleRate.saveSettingsTo(settings);
		m_balanceClassInstances.saveSettingsTo(settings);

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_imgColumn.validateSettings(settings);
		m_labelingColumn.validateSettings(settings);
		m_dimSelection.validateSettings(settings);
		m_featDimSelection.validateSettings(settings);
		m_classifierSelection.validateSettings(settings);
		m_resampleRate.validateSettings(settings);
		m_balanceClassInstances.validateSettings(settings);

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_imgColumn.loadSettingsFrom(settings);
		m_labelingColumn.loadSettingsFrom(settings);
		m_dimSelection.loadSettingsFrom(settings);
		m_featDimSelection.loadSettingsFrom(settings);
		m_classifierSelection.loadSettingsFrom(settings);
		m_resampleRate.loadSettingsFrom(settings);
		m_balanceClassInstances.loadSettingsFrom(settings);

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		//

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec)
			throws IOException, CanceledExecutionException {
		//

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec)
			throws IOException, CanceledExecutionException {
		//

	}

}
