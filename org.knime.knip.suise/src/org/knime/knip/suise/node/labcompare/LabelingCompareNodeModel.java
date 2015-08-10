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
package org.knime.knip.suise.node.labcompare;

import java.io.File;
import java.io.IOException;

import net.imagej.ImgPlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.numeric.integer.ByteType;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.RowKey;
import org.knime.core.data.container.CellFactory;
import org.knime.core.data.container.ColumnRearranger;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.knip.base.data.img.ImgPlusCell;
import org.knime.knip.base.data.img.ImgPlusCellFactory;
import org.knime.knip.base.data.labeling.LabelingValue;
import org.knime.knip.base.node.NodeUtils;
import org.knime.knip.core.KNIPGateway;
import org.knime.knip.core.features.FeatureFactory;
import org.knime.knip.core.features.FeatureSet;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class LabelingCompareNodeModel<L extends Comparable<L>> extends NodeModel {

	private SettingsModelString m_labelingCol1 = createLabelingCol1Model();

	private SettingsModelString m_labelingCol2 = createLabelingCol2Model();

	static SettingsModelString createLabelingCol1Model() {
		return new SettingsModelString("labeling_column_1", "");
	}

	static SettingsModelString createLabelingCol2Model() {
		return new SettingsModelString("labeling_column_2", "");
	}

	/**
	 * @param nrInDataPorts
	 * @param nrOutDataPorts
	 */
	protected LabelingCompareNodeModel() {
		super(1, 1);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(DataTableSpec[] inSpecs) throws InvalidSettingsException {
		NodeUtils.autoColumnSelection(inSpecs[0], m_labelingCol1, LabelingValue.class, this.getClass());
		NodeUtils.autoColumnSelection(inSpecs[0], m_labelingCol2, LabelingValue.class, this.getClass());
		ColumnRearranger rearranger = new ColumnRearranger(inSpecs[0]);
		rearranger.append(createCellFactory(0, 1, null, new LabelingCompareFeatureSet<L>()));
		return new DataTableSpec[] { rearranger.createSpec() };

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected BufferedDataTable[] execute(BufferedDataTable[] inData, ExecutionContext exec) throws Exception {

		int colIdx1 = NodeUtils.autoColumnSelection(inData[0].getDataTableSpec(), m_labelingCol1, LabelingValue.class,
				this.getClass());
		int colIdx2 = NodeUtils.autoColumnSelection(inData[0].getDataTableSpec(), m_labelingCol2, LabelingValue.class,
				this.getClass());

		ColumnRearranger rearranger = new ColumnRearranger(inData[0].getDataTableSpec());
		rearranger.append(createCellFactory(colIdx1, colIdx2, exec, new LabelingCompareFeatureSet<L>()));

		return new BufferedDataTable[] { exec.createColumnRearrangeTable(inData[0], rearranger, exec) };

	}

	private CellFactory createCellFactory(final int colIdx1, final int colIdx2, ExecutionContext exec,
			final FeatureSet... featureSets) {

		final FeatureFactory featFac = new FeatureFactory(true, featureSets);

		final ImgPlusCellFactory cellFac = exec != null ? new ImgPlusCellFactory(exec) : null;

		return new CellFactory() {

			@Override
			public void setProgress(int curRowNr, int rowCount, RowKey lastKey, ExecutionMonitor exec) {
				exec.setProgress(curRowNr / (double) rowCount);

			}

			@Override
			public DataColumnSpec[] getColumnSpecs() {
				DataColumnSpec[] specs = new DataColumnSpec[featFac.getNumFeatures() + 1];

				specs[0] = new DataColumnSpecCreator("Difference Image", ImgPlusCell.TYPE).createSpec();
				for (int i = 1; i < specs.length; i++) {
					specs[i] = new DataColumnSpecCreator(featFac.getFeatureNames()[i - 1], DoubleCell.TYPE)
							.createSpec();
				}

				return specs;
			}

			@Override
			public DataCell[] getCells(DataRow row) {
				if (row.getCell(colIdx1).isMissing()) {
					throw new RuntimeException("Reference labeling is missing (missing cell)!");
				}

				RandomAccessibleInterval<LabelingType<L>>[] labs = new RandomAccessibleInterval[2];
				labs[0] = ((LabelingValue<L>) row.getCell(colIdx1)).getLabeling();
				if (row.getCell(colIdx2).isMissing()) {
					// create empty labeling
					setWarningMessage("Missing target labeling (missing cell). Empty labeling assumed.");
					labs[1] = KNIPGateway.ops().create().imgLabeling(labs[0]);
				} else {
					labs[1] = ((LabelingValue<L>) row.getCell(colIdx2)).getLabeling();
				}
				featFac.updateFeatureTarget(labs);

				DataCell[] cells = new DataCell[featFac.getNumFeatures() + 1];
				for (FeatureSet fset : featureSets) {
					if (fset instanceof LabelingCompareFeatureSet) {
						Img<ByteType> img = ((LabelingCompareFeatureSet) fset).getDiffImg();
						try {
							cells[0] = cellFac.createCell(new ImgPlus<ByteType>(img));
						} catch (IOException e) {
							throw new RuntimeException(e);
						}
					}
				}
				if (cells[0] == null) {
					cells[0] = DataType.getMissingCell();
				}

				for (int i = 1; i < cells.length; i++) {
					double val = featFac.getFeatureValue(i - 1);
					if (Double.isNaN(val)) {
						cells[i] = DataType.getMissingCell();
					} else {
						cells[i] = new DoubleCell(val);
					}
				}
				return cells;
			}
		};
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec)
			throws IOException, CanceledExecutionException {
		// TODO Auto-generated method stub

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec)
			throws IOException, CanceledExecutionException {
		// TODO Auto-generated method stub

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(NodeSettingsWO settings) {
		m_labelingCol1.saveSettingsTo(settings);
		m_labelingCol2.saveSettingsTo(settings);

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		m_labelingCol1.validateSettings(settings);
		m_labelingCol1.validateSettings(settings);

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
		m_labelingCol1.loadSettingsFrom(settings);
		m_labelingCol2.loadSettingsFrom(settings);

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		// TODO Auto-generated method stub

	}

}
