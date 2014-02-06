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

import java.util.List;

import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.meta.ImgPlus;
import net.imglib2.ops.operation.randomaccessibleinterval.unary.regiongrowing.AbstractRegionGrowing;
import net.imglib2.ops.types.ConnectedType;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.integer.IntType;

import org.knime.core.node.ExecutionContext;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.data.labeling.LabelingCell;
import org.knime.knip.base.data.labeling.LabelingCellFactory;
import org.knime.knip.base.node.ValueToCellNodeDialog;
import org.knime.knip.base.node.ValueToCellNodeFactory;
import org.knime.knip.base.node.ValueToCellNodeModel;
import org.knime.knip.core.data.img.DefaultLabelingMetadata;
import org.knime.knip.core.util.EnumUtils;

/**
 * Thresholds the input image at those intensity levels, which result in
 * different binary images and performs a connected component analysis on each.
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class LevelCCANodeFactory<T extends IntegerType<T> & NativeType<T>>
        extends ValueToCellNodeFactory<ImgPlusValue<T>> {

    private static SettingsModelString createTypeModel() {
        return new SettingsModelString("connection_type",
                ConnectedType.values()[0].toString());
    }

    private static SettingsModelBoolean createBackgroundModel() {
        return new SettingsModelBoolean("white_background", false);
    }

    private static SettingsModelIntegerBounded createLowerBoundModel() {
        return new SettingsModelIntegerBounded("lower_bound", 0, 0, 255);
    }

    private static SettingsModelIntegerBounded createUpperBoundModel() {
        return new SettingsModelIntegerBounded("upper_bound", 255, 0, 255);
    }

    private static SettingsModelIntegerBounded createStepSizeModel() {
        return new SettingsModelIntegerBounded("step_size", 1, 0, 255);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ValueToCellNodeModel<ImgPlusValue<T>, LabelingCell<Integer>> createNodeModel() {
        return new ValueToCellNodeModel<ImgPlusValue<T>, LabelingCell<Integer>>() {

            private LabelingCellFactory m_labelingCellFactory;

            private SettingsModelString m_type = createTypeModel();

            private SettingsModelBoolean m_background = createBackgroundModel();

            private SettingsModelInteger m_lowerBound = createLowerBoundModel();

            private SettingsModelInteger m_upperBound = createUpperBoundModel();

            private SettingsModelIntegerBounded m_stepSize =
                    createStepSizeModel();

            /**
             * {@inheritDoc}
             */
            @Override
            protected void prepareExecute(ExecutionContext exec) {
                m_labelingCellFactory = new LabelingCellFactory(exec);
            }

            @Override
            protected LabelingCell<Integer> compute(ImgPlusValue<T> cellValue)
                    throws Exception {
                ImgPlus<T> img = cellValue.getImgPlus();

                long[][] structuringElement;
                if (m_type.getStringValue().equals(
                        ConnectedType.EIGHT_CONNECTED.name())) {
                    structuringElement =
                            AbstractRegionGrowing.get8ConStructuringElement(img
                                    .numDimensions());
                } else {
                    structuringElement =
                            AbstractRegionGrowing.get4ConStructuringElement(img
                                    .numDimensions());
                }

                Labeling<Integer> res =
                        new NativeImgLabeling<Integer, IntType>(
                                new ArrayImgFactory<IntType>().create(img,
                                        new IntType()));
                new LevelCCA<T>(structuringElement,
                        m_background.getBooleanValue(),
                        m_lowerBound.getIntValue(), m_upperBound.getIntValue(),
                        m_stepSize.getIntValue()).compute(img, res);
                return m_labelingCellFactory.createCell(res,
                        new DefaultLabelingMetadata(img, img, img, null));
            }

            @Override
            protected void addSettingsModels(List<SettingsModel> settingsModels) {
                settingsModels.add(m_background);
                settingsModels.add(m_type);
                settingsModels.add(m_lowerBound);
                settingsModels.add(m_upperBound);
                settingsModels.add(m_stepSize);
            }
        };
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected ValueToCellNodeDialog<ImgPlusValue<T>> createNodeDialog() {
        return new ValueToCellNodeDialog<ImgPlusValue<T>>() {

            @Override
            public void addDialogComponents() {
                addDialogComponent(
                        "Options",
                        "Settings",
                        new DialogComponentStringSelection(
                                createTypeModel(),
                                "Connection Type",
                                EnumUtils
                                        .getStringListFromToString(ConnectedType
                                                .values())));
                addDialogComponent("Options", "Settings",
                        new DialogComponentBoolean(createBackgroundModel(),
                                "White background"));

                addDialogComponent(
                        "Options",
                        "Settings",
                        new DialogComponentNumber(createLowerBoundModel(),
                                "Lower pixel value bound (0-255, inclusive)", 1));

                addDialogComponent(
                        "Options",
                        "Settings",
                        new DialogComponentNumber(createUpperBoundModel(),
                                "Upper pixel value bound (0-255, inclusive)", 1));

                addDialogComponent("Options", "Settings",
                        new DialogComponentNumber(createStepSizeModel(),
                                "Step size", 1));
            }
        };
    }
}
