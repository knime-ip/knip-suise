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

import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.port.PortObjectSpec;

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Classifier;

/**
 * 
 * @author hornm, University of Konstanz
 */
public class SettingsModelWekaClassifier extends SettingsModel {

    private static final String KEY_CLASSIFIER_NAME = "_classifier_name";

    private static final String KEY_CLASSIFIER_OPTIONS = "_classifier_options";

    private final String m_configName;

    private AbstractClassifier m_wekaClassifier;

    /**
     * Creates a new model.
     * 
     * @param configName
     * @param defaultClassifier
     */
    public SettingsModelWekaClassifier(String configName,
            AbstractClassifier defaultClassifier) {
        m_configName = configName;
        m_wekaClassifier = defaultClassifier;

    }

    /**
     * @return the classifier
     */
    public Classifier getClassifier() {
        return m_wekaClassifier;
    }

    /**
     * @param classifier the new classifier
     */
    public void setClassifier(AbstractClassifier classifier) {
        m_wekaClassifier = classifier;
    }

    /**
     * {@inheritDoc}
     */
    @SuppressWarnings("unchecked")
    @Override
    protected SettingsModelWekaClassifier createClone() {
        return new SettingsModelWekaClassifier(m_configName, m_wekaClassifier);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String getModelTypeID() {
        return "SMID_weka_classifier";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String getConfigName() {
        return m_configName;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadSettingsForDialog(NodeSettingsRO settings,
            PortObjectSpec[] specs) throws NotConfigurableException {
        try {
            loadSettingsForModel(settings);
        } catch (InvalidSettingsException e) {
            // keep default classifier
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsForDialog(NodeSettingsWO settings)
            throws InvalidSettingsException {
        saveSettingsForModel(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettingsForModel(NodeSettingsRO settings)
            throws InvalidSettingsException {

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadSettingsForModel(NodeSettingsRO settings)
            throws InvalidSettingsException {
        String name = settings.getString(m_configName + KEY_CLASSIFIER_NAME);
        String[] options =
                settings.getStringArray(m_configName + KEY_CLASSIFIER_OPTIONS);
        try {
            m_wekaClassifier =
                    (AbstractClassifier)AbstractClassifier.forName(name,
                            options);
        } catch (Exception e) {
            throw new InvalidSettingsException(e);
        }

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsForModel(NodeSettingsWO settings) {
        settings.addString(m_configName + KEY_CLASSIFIER_NAME, m_wekaClassifier
                .getClass().getCanonicalName());
        settings.addStringArray(m_configName + KEY_CLASSIFIER_OPTIONS,
                m_wekaClassifier.getOptions());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return m_wekaClassifier.toString();
    }

}
