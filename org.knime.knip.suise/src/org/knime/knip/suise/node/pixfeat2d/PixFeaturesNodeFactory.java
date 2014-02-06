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
import java.util.List;

import net.imglib2.type.numeric.RealType;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;
import org.knime.core.node.defaultnodesettings.DialogComponent;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.DialogComponentLabel;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.node.ValueToCellNodeModel;
import org.knime.knip.base.nodes.features.FeatureSetDialogComponentCollection;
import org.knime.knip.base.nodes.view.TableCellViewNodeView;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public abstract class PixFeaturesNodeFactory<T extends RealType<T>> extends
        NodeFactory<PixFeaturesNodeModel<T>> {

    /**
     * @return
     */
    protected abstract PixFeatureSetProvider[] getPixFeatureSetProviders();

    protected abstract int getPixFeatureDimensionality();

    protected abstract int getMaxNumOrientations();

    /**
     * {@inheritDoc}
     */
    @Override
    protected int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected NodeDialogPane createNodeDialogPane() {
        FeatureSetDialogComponentCollection featColl =
                new FeatureSetDialogComponentCollection(
                        PixFeaturesNodeModel.createActiveFeatureSetModel());

        List<DialogComponent> diagComps = new ArrayList<DialogComponent>();

        // feature dialog components
        for (PixFeatureSetProvider<T> p : getPixFeatureSetProviders()) {
            p.initAndAddDialogComponents(diagComps);
            featColl.addFeatureSetDialogComponent(p.getFeatureSetId(),
                    p.getFeatureSetName(),
                    new DialogComponentLabel(p.getFeatureSetName()));
            for (DialogComponent dc : diagComps) {
                featColl.addFeatureSetDialogComponent(p.getFeatureSetId(),
                        p.getFeatureSetName(), dc);
            }
            diagComps.clear();
        }

        // dialog components for feature selection
        featColl.addDialogComponent(
                "Column Selection",
                "Creation Mode",
                new DialogComponentStringSelection(ValueToCellNodeModel
                        .createColCreationModeModel(), "Column Creation Mode",
                        ValueToCellNodeModel.COL_CREATION_MODES));
        featColl.addDialogComponent(
                "Column Selection",
                "Column suffix",
                new DialogComponentString(ValueToCellNodeModel
                        .createColSuffixNodeModel(), "Column suffix"));

        featColl.addDialogComponent(
                "Column Selection",
                "",
                new DialogComponentColumnFilter(ValueToCellNodeModel
                        .createColumnSelectionModel(), 0, false,
                        ImgPlusValue.class));

        SettingsModelIntegerBounded maxNumOrientations =
                PixFeaturesNodeModel
                        .createNumOrientationsModel(getMaxNumOrientations());
        SettingsModelString orientationDimLabel =
                PixFeaturesNodeModel.createOrientationDimLabelModel();
        maxNumOrientations.setIntValue(getMaxNumOrientations());
        if (maxNumOrientations.getIntValue() == 1) {
            maxNumOrientations.setEnabled(false);
            orientationDimLabel.setEnabled(false);
        }
        featColl.addDialogComponent("Additional settings", "",
                new DialogComponentNumber(maxNumOrientations,
                        "Number of orientations", 1));
        featColl.addDialogComponent("Additional settings", "",
                new DialogComponentString(orientationDimLabel,
                        "Orientation dimension label"));
        featColl.addDialogComponent(
                "Additional settings",
                "",
                new DialogComponentString(PixFeaturesNodeModel
                        .createFeatureDimLabelModel(),
                        "Feature dimension label"));

        return featColl.getDialog();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public PixFeaturesNodeModel<T> createNodeModel() {
        return new PixFeaturesNodeModel<T>(getPixFeatureSetProviders(),
                getPixFeatureDimensionality(), getMaxNumOrientations());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<PixFeaturesNodeModel<T>> createNodeView(int viewIndex,
            PixFeaturesNodeModel<T> nodeModel) {
        return new TableCellViewNodeView<PixFeaturesNodeModel<T>>(nodeModel);
    }

}
