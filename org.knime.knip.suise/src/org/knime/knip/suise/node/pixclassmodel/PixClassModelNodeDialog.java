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

import java.util.Set;

import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnNameSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.data.labeling.LabelingValue;
import org.knime.knip.base.node.dialog.DialogComponentDimSelection;
import org.knime.knip.base.node.nodesettings.SettingsModelDimSelection;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class PixClassModelNodeDialog extends DefaultNodeSettingsPane {

    private SettingsModelDimSelection m_smImageDims;

    private SettingsModelDimSelection m_smFeatureDim;

    public PixClassModelNodeDialog() {

        createNewGroup("Classifier selection");
        addDialogComponent(new DialogComponentWekaClassifier(
                PixClassModelNodeModel.createClassifierSelectionModel(), ""));
        closeCurrentGroup();

        createNewGroup("Column selection");
        addDialogComponent(new DialogComponentColumnNameSelection(
                PixClassModelNodeModel.createImgColumnModel(), "Feature image",
                PixClassModelNodeModel.INPUT_DATATABLE_PORT, ImgPlusValue.class));

        addDialogComponent(new DialogComponentColumnNameSelection(
                PixClassModelNodeModel.createLabColumnModel(), "Labeling",
                PixClassModelNodeModel.INPUT_DATATABLE_PORT,
                LabelingValue.class));
        closeCurrentGroup();

        createNewGroup("Dimension selection");
        m_smImageDims = PixClassModelNodeModel.createDimSelectionModel();
        addDialogComponent(new DialogComponentDimSelection(m_smImageDims,
                "Image dimensions", 2, 5));

        m_smFeatureDim = PixClassModelNodeModel.createFeatDimSelectionModel();
        addDialogComponent(new DialogComponentDimSelection(m_smFeatureDim,
                "Feature dimension", 1, 1));
        closeCurrentGroup();

        createNewGroup("Additional options");
        addDialogComponent(new DialogComponentNumber(
                PixClassModelNodeModel.createSampleRateModel(),
                "Resample Rate", .01));

        addDialogComponent(new DialogComponentBoolean(
                PixClassModelNodeModel.createBalanceModel(),
                "Balance training set"));
        closeCurrentGroup();

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void saveAdditionalSettingsTo(NodeSettingsWO settings)
            throws InvalidSettingsException {
        super.saveAdditionalSettingsTo(settings);

        // check that the feature dimension is not already selected as image
        // dimension
        Set<String> featLabels = m_smFeatureDim.getSelectedDimLabels();
        Set<String> imageLabels = m_smImageDims.getSelectedDimLabels();
        if (imageLabels.contains(featLabels.iterator().next())) {
            throw new InvalidSettingsException(
                    "Selected feature dimensions must not be an image dimension!");
        }

    }
}
