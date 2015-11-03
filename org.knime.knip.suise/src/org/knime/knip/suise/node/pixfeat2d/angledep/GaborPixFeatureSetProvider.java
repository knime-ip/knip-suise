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
package org.knime.knip.suise.node.pixfeat2d.angledep;

import java.util.ArrayList;
import java.util.List;

import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.DoubleType;

import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.DialogComponent;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.knip.suise.data.feat.PixFeatureSet;
import org.knime.knip.suise.node.pixfeat2d.PixFeatureSetProvider;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class GaborPixFeatureSetProvider<T extends RealType<T>> extends
        PixFeatureSetProvider<T> {

    private static NodeLogger LOGGER = NodeLogger
            .getLogger(GaborPixFeatureSetProvider.class);

    private static final SettingsModelInteger createRadiusModel() {
        return new SettingsModelInteger("radius", 30);
    }

    private static final SettingsModelString createScalesModel() {
        return new SettingsModelString("scales", "0.5,1,1.5");
    }

    private static final SettingsModelString createFrequenciesModel() {
        return new SettingsModelString("frequencies", "1");
    }

    private static final SettingsModelString createElongationsModel() {
        return new SettingsModelString("elongations", "2");
    }

    private SettingsModelInteger m_smRadius;

    private SettingsModelString m_smScales;

    private SettingsModelString m_smFrequencies;

    private SettingsModelString m_smElongations;

    /**
     * {@inheritDoc}
     */
    @Override
    protected PixFeatureSet<T> getPixFeatureSet(int numOrientations) {
        return new GaborFilterFeatureSet<T, DoubleType>(
                extractDoubles(m_smScales.getStringValue()),
                extractDoubles(m_smFrequencies.getStringValue()),
                extractDoubles(m_smElongations.getStringValue()),
                m_smRadius.getIntValue(), numOrientations, false,
                new DoubleType());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void initAndAddSettingsModels(List<SettingsModel> settingsModels) {
        settingsModels.add(m_smRadius = createRadiusModel());
        settingsModels.add(m_smScales = createScalesModel());
        settingsModels.add(m_smFrequencies = createFrequenciesModel());
        settingsModels.add(m_smElongations = createElongationsModel());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void initAndAddDialogComponents(
            List<DialogComponent> dialogComponents) {
        DialogComponent dc =
                new DialogComponentNumber(createRadiusModel(), "Radius", 1);
        dialogComponents.add(dc);
        dc = new DialogComponentString(createScalesModel(), "Scales");
        dialogComponents.add(dc);
        dc = new DialogComponentString(createFrequenciesModel(), "Frequencies");
        dialogComponents.add(dc);
        dc = new DialogComponentString(createElongationsModel(), "Elongations");
        dialogComponents.add(dc);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getFeatureSetName() {
        return "Gabor feature set";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getFeatureSetId() {
        return "Gabor feature set";
    }

    private double[] extractDoubles(String s) {
        String[] split = s.split(",");
        ArrayList<Double> tmp = new ArrayList<Double>(split.length);
        for (String n : split) {
            try {
                tmp.add(Double.parseDouble(n));
            } catch (NumberFormatException e) {
                LOGGER.warn("Cannot parse number from " + n);
            }
        }
        double[] res = new double[tmp.size()];
        for (int i = 0; i < res.length; i++) {
            res[i] = tmp.get(i);
        }
        return res;

    }

	@Override
	public void cleanUp() {
		// empty
	}

}
