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
package org.knime.knip.suise.node.pixfeat2d.fiji;

import ij.ImagePlus;
import ij.process.ImageProcessor;

import java.util.List;

import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.meta.Axes;
import net.imglib2.meta.CalibratedAxis;
import net.imglib2.meta.ImgPlus;
import net.imglib2.meta.ImgPlusMetadata;
import net.imglib2.meta.axis.DefaultLinearAxis;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

import org.knime.core.node.ExecutionContext;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;
import org.knime.knip.base.data.img.ImgPlusCell;
import org.knime.knip.base.data.img.ImgPlusCellFactory;
import org.knime.knip.base.data.img.ImgPlusValue;
import org.knime.knip.base.node.ValueToCellNodeModel;

import trainableSegmentation.FeatureStack;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class FijiTrainSegFeaturesNodeModel<T extends RealType<T>> extends
        ValueToCellNodeModel<ImgPlusValue<T>, ImgPlusCell<FloatType>> {

    final static SettingsModelStringArray createFeatureListModel() {
        // default selection
        String[] def = new String[5];
        for (int i = 0; i < def.length; i++) {
            def[i] = FeatureStack.availableFeatures[i];
        }
        return new SettingsModelStringArray("feature_list", def);
    }

    final static SettingsModelInteger createMembraneThicknessModel() {
        return new SettingsModelInteger("membran_thickness", 1);
    }

    final static SettingsModelInteger createMembranePatchSizeModel() {
        return new SettingsModelInteger("membrane_patch_size", 19);
    }

    final static SettingsModelInteger createMinSigmaModel() {
        return new SettingsModelInteger("min_sigma", 1);
    }

    final static SettingsModelInteger createMaxSigmaModel() {
        return new SettingsModelInteger("max_sigma", 16);
    }

    final static SettingsModelBoolean createMultiThreadedModel() {
        return new SettingsModelBoolean("multi_threaded", true);
    }

    final static SettingsModelString createFeatDimLabelModel() {
        return new SettingsModelString("feature_dim_label", "F");
    }

    private SettingsModelStringArray m_smFeatureList = createFeatureListModel();

    private SettingsModelInteger m_smMembraneThickness =
            createMembraneThicknessModel();

    private SettingsModelInteger m_smMembranePatchSize =
            createMembranePatchSizeModel();

    private SettingsModelInteger m_smMinSigma = createMinSigmaModel();

    private SettingsModelInteger m_smMaxSigma = createMaxSigmaModel();

    private SettingsModelBoolean m_smMultiThreaded = createMultiThreadedModel();

    private SettingsModelString m_smFeatDimLabel = createFeatDimLabelModel();

    private boolean[] m_enabledFeatures;

    private ImgPlusCellFactory m_imgCellFactory;

    /**
     * {@inheritDoc}
     */
    @Override
    protected void addSettingsModels(List<SettingsModel> settingsModels) {
        settingsModels.add(m_smFeatureList);
        settingsModels.add(m_smMembraneThickness);
        settingsModels.add(m_smMembranePatchSize);
        settingsModels.add(m_smMinSigma);
        settingsModels.add(m_smMaxSigma);
        settingsModels.add(m_smMultiThreaded);
        settingsModels.add(m_smFeatDimLabel);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void prepareExecute(ExecutionContext exec) {
        m_enabledFeatures = new boolean[FeatureStack.availableFeatures.length];
        String[] selected = m_smFeatureList.getStringArrayValue();
        int j = 0;
        for (int i = 0; i < FeatureStack.availableFeatures.length
                && j < selected.length; i++) {
            if (FeatureStack.availableFeatures[i].equals(selected[j])) {
                m_enabledFeatures[i] = true;
                j++;
            }
        }
        m_imgCellFactory = new ImgPlusCellFactory(exec);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected ImgPlusCell<FloatType> compute(ImgPlusValue<T> cellValue)
            throws Exception {
        ImgPlus<T> img = cellValue.getImgPlus();

        // create feature stack
        if (img.numDimensions() > 2) {
            throw new IllegalArgumentException("Only 2D images supported, yet!");
        }

        ImagePlus ip = new ImagePlus();
        new ImgToIJ(2).compute(img, ip);

        FeatureStack ijFeatureStack = new FeatureStack(ip);

        // set parameters
        ijFeatureStack.setEnabledFeatures(m_enabledFeatures);
        ijFeatureStack
                .setMembranePatchSize(m_smMembranePatchSize.getIntValue());
        ijFeatureStack.setMembraneSize(m_smMembraneThickness.getIntValue());
        ijFeatureStack.setMinimumSigma(m_smMinSigma.getIntValue());
        ijFeatureStack.setMaximumSigma(m_smMaxSigma.getIntValue());

        // caluclate all features
        if (m_smMultiThreaded.getBooleanValue()) {
            ijFeatureStack.updateFeaturesMT();
        } else {
            ijFeatureStack.updateFeaturesST();
        }

        // convert image stack to img
        Img<FloatType> res =
                new ArrayImgFactory<FloatType>().create(
                        new long[]{img.dimension(0), img.dimension(1),
                                ijFeatureStack.getSize()}, new FloatType());
        Cursor<FloatType> resCur = res.cursor();
        for (int i = 0; i < ijFeatureStack.getSize(); i++) {
            for (int j = 0; j < img.size(); j++) {
                resCur.fwd();
                ImageProcessor tmp = ijFeatureStack.getProcessor(i + 1);
                resCur.get().set(tmp.get(j));
            }
        }

        CalibratedAxis[] axes =
                new CalibratedAxis[]{
                        new DefaultLinearAxis(Axes.X),
                        new DefaultLinearAxis(Axes.Y),
                        new DefaultLinearAxis(Axes.get(m_smFeatDimLabel
                                .getStringValue()))};
        ImgPlusMetadata metadata = cellValue.getMetadata();
        for (int i = 0; i < axes.length; i++) {
            metadata.setAxis(axes[i], i);
        }

        return m_imgCellFactory
                .createCell(new ImgPlus<FloatType>(res, metadata));
    }
}
