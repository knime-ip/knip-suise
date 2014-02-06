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
package org.knime.knip.suise.node.pixclassapply;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.imglib2.img.Img;
import net.imglib2.meta.Axes;
import net.imglib2.meta.CalibratedAxis;
import net.imglib2.meta.DefaultCalibratedSpace;
import net.imglib2.meta.ImgPlus;
import net.imglib2.meta.axis.DefaultLinearAxis;
import net.imglib2.ops.operation.Operations;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.ShortType;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowIterator;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.BufferedDataTableHolder;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.PortObjectSpec;
import org.knime.core.node.port.PortType;
import org.knime.knip.base.data.img.ImgPlusCell;
import org.knime.knip.base.data.img.ImgPlusCellFactory;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.data.labeling.LabelingCell;
import org.knime.knip.base.data.labeling.LabelingCellFactory;
import org.knime.knip.base.node.NodeUtils;
import org.knime.knip.base.node.nodesettings.SettingsModelDimSelection;
import org.knime.knip.core.data.img.DefaultImgMetadata;
import org.knime.knip.core.data.img.DefaultLabelingMetadata;
import org.knime.knip.suise.node.port.WekaClassifierPortObject;
import org.knime.knip.suise.node.port.WekaClassifierPortObjectSpec;
import org.knime.knip.suise.ops.ClassifyPixels;

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.Classifier;
import weka.core.Instances;

/**
 * 
 * @author hornm, University of Konstanz
 */
public class PixClassApplyNodeModel<T extends RealType<T>, L extends Comparable<L>>
        extends NodeModel implements BufferedDataTableHolder {

    private static final NodeLogger LOGGER = NodeLogger
            .getLogger(PixClassApplyNodeModel.class);

    static final SettingsModelString createResultSelectionModel() {
        return new SettingsModelString("result_type_selection",
                RESULT_CALCULATION[0]);
    }

    static final String[] RESULT_CALCULATION = new String[]{
            "Probability Maps AND Labeling", "Probability Maps", "Labeling"};

    static final SettingsModelString createImgColumnModel() {
        return new SettingsModelString("img_column", "");
    }

    static final SettingsModelDimSelection createDimSelectionModel() {
        return new SettingsModelDimSelection("dim_selection", "X", "Y");
    }

    static final SettingsModelDimSelection createFeatDimSelectionModel() {
        return new SettingsModelDimSelection("feat_dim_selection", "C");
    }

    static final int INPUT_DATATABLE_PORT = 1;

    private SettingsModelString m_imgColumn = createImgColumnModel();

    private SettingsModelDimSelection m_dimSelection =
            createDimSelectionModel();

    private SettingsModelDimSelection m_featDimSelection =
            createFeatDimSelectionModel();

    private SettingsModelString m_resultCalc = createResultSelectionModel();

    // data for the table cell view
    private BufferedDataTable m_data;

    protected PixClassApplyNodeModel() {
        super(new PortType[]{WekaClassifierPortObject.TYPE,
                BufferedDataTable.TYPE}, new PortType[]{BufferedDataTable.TYPE});
    }

    /**
     * {@inheritDoc}
     */
    protected PortObjectSpec[] configure(PortObjectSpec[] inSpecs)
            throws InvalidSettingsException {
        getImgColumnIndex((DataTableSpec)inSpecs[INPUT_DATATABLE_PORT]);
        return new PortObjectSpec[]{createOutSpec((WekaClassifierPortObjectSpec)inSpecs[0])};
    }

    /**
     * {@inheritDoc}
     */
    protected PortObject[] execute(PortObject[] inObjects, ExecutionContext exec)
            throws Exception {

        BufferedDataTable inTable =
                (BufferedDataTable)inObjects[INPUT_DATATABLE_PORT];

        int imgColIdx = getImgColumnIndex(inTable.getDataTableSpec());

        Classifier classifier =
                ((WekaClassifierPortObject)inObjects[0]).getClassifier();
        Instances trainingSet =
                ((WekaClassifierPortObject)inObjects[0]).getTrainingInstances();

        List<String> labels =
                Arrays.asList(((WekaClassifierPortObject)inObjects[0])
                        .getSpec().getClassLabels());

        /* classification process */
        BufferedDataContainer resContainer =
                exec.createDataContainer(createOutSpec(((WekaClassifierPortObject)inObjects[0])
                        .getSpec()));

        RowIterator it;
        DataRow row;
        it = inTable.iterator();
        int rowCount = inTable.getRowCount();
        int i = 0;
        ImgPlusCellFactory imgCellFactory = new ImgPlusCellFactory(exec);
        LabelingCellFactory labCellFactory = new LabelingCellFactory(exec);
        List<DataCell> resCells = new ArrayList<DataCell>();
        while (it.hasNext()) {
            row = it.next();
            ImgPlus<T> img =
                    ((ImgPlusValue<T>)row.getCell(imgColIdx)).getImgPlus();
            int[] tmp =
                    m_featDimSelection.getSelectedDimIndices(
                            img.numDimensions(), img);
            if (tmp.length == 0) {
                LOGGER.warn("Feature dimensions doesn't exist in image in row "
                        + row.getKey() + ".");
                continue;
            }
            if (img.dimension(tmp[0]) != trainingSet.numAttributes() - 1) {
                LOGGER.warn("The size of the feature dimensions for image in row "
                        + row.getKey()
                        + " doesn't equal the number of features of the classifer. Row ignored!");
                continue;
            }
            int featDim = tmp[0];

            int[] dimIndices =
                    m_dimSelection.getSelectedDimIndices(img.numDimensions(),
                            img);

            // classify image
            exec.setProgress("Classifing image in row " + row.getKey());
            ClassifyPixels<T, ShortType> cp =
                    new ClassifyPixels<T, ShortType>(classifier, trainingSet,
                            labels.size(), dimIndices, featDim, !m_resultCalc
                                    .getStringValue().equals(
                                            RESULT_CALCULATION[1]),
                            new ShortType());
            Img<ShortType>[] probs = Operations.compute(cp, img);

            // create results (probabilities and/or labelings)
            CalibratedAxis[] axes = new CalibratedAxis[dimIndices.length];
            for (int j = 0; j < axes.length; j++) {
                axes[j] =
                        new DefaultLinearAxis(Axes.get(img.axis(dimIndices[j])
                                .type().getLabel()));
            }
            DefaultImgMetadata metadata =
                    new DefaultImgMetadata(new DefaultCalibratedSpace(axes),
                            img, img, img);
            resCells.clear();
            if (!m_resultCalc.getStringValue().equals(RESULT_CALCULATION[2])) {
                for (int j = 0; j < probs.length; j++) {
                    resCells.add(imgCellFactory
                            .createCell(new ImgPlus<ShortType>(probs[j],
                                    metadata)));
                }
            }
            if (!m_resultCalc.getStringValue().equals(RESULT_CALCULATION[1])) {
                resCells.add(labCellFactory.createCell(cp.getLabeling(),
                        new DefaultLabelingMetadata(img, img, img, null)));
            }
            resContainer.addRowToTable(new DefaultRow(row.getKey(), resCells));
            exec.checkCanceled();
            exec.setProgress((double)i / rowCount);
            i++;
        }
        resContainer.close();
        m_data = resContainer.getTable();
        return new BufferedDataTable[]{m_data};

    }

    private final int getImgColumnIndex(DataTableSpec inSpec)
            throws InvalidSettingsException {
        return NodeUtils.silentOptionalAutoColumnSelection(inSpec, m_imgColumn,
                ImgPlusValue.class);

    }

    private final DataTableSpec createOutSpec(
            WekaClassifierPortObjectSpec wcSpec) {
        List<DataColumnSpec> resSpec = new ArrayList<DataColumnSpec>();
        String[] labels = wcSpec.getClassLabels();
        if (!m_resultCalc.getStringValue().equals(RESULT_CALCULATION[2])) {
            for (String label : labels) {
                resSpec.add(new DataColumnSpecCreator(label, ImgPlusCell.TYPE)
                        .createSpec());
            }
        }
        if (!m_resultCalc.getStringValue().equals(RESULT_CALCULATION[1])) {
            resSpec.add(new DataColumnSpecCreator("Predicted labeling",
                    LabelingCell.TYPE).createSpec());
        }
        return new DataTableSpec(resSpec.toArray(new DataColumnSpec[resSpec
                .size()]));

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(NodeSettingsWO settings) {
        m_imgColumn.saveSettingsTo(settings);
        m_dimSelection.saveSettingsTo(settings);
        m_featDimSelection.saveSettingsTo(settings);
        m_resultCalc.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_imgColumn.validateSettings(settings);
        m_dimSelection.validateSettings(settings);
        m_featDimSelection.validateSettings(settings);
        m_resultCalc.validateSettings(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_imgColumn.loadSettingsFrom(settings);
        m_dimSelection.loadSettingsFrom(settings);
        m_featDimSelection.loadSettingsFrom(settings);
        m_resultCalc.loadSettingsFrom(settings);

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

    // Necessary for the Table cell view

    /**
     * {@inheritDoc}
     */
    @Override
    public BufferedDataTable[] getInternalTables() {
        return new BufferedDataTable[]{m_data};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setInternalTables(BufferedDataTable[] tables) {
        m_data = tables[0];

    }

}
