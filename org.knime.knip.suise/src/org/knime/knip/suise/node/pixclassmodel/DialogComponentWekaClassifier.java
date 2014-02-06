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

import java.awt.FlowLayout;
import java.util.Properties;

import javax.swing.JLabel;
import javax.swing.JPanel;

import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DialogComponent;
import org.knime.core.node.port.PortObjectSpec;

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Classifier;
import weka.core.Environment;
import weka.core.Utils;
import weka.gui.GenericObjectEditor;
import weka.gui.PropertyPanel;
import weka.gui.beans.PluginManager;

/**
 * 
 * @author Martin Horn, University of Konstanz
 */
public class DialogComponentWekaClassifier extends DialogComponent {

    static {
        // deactivate the package management of weka
        Environment env = Environment.getSystemWide();
        env.addVariableSystemWide("weka.core.loadPackages", "false");
        env.addVariableSystemWide("weka.packageManager.offline", "true");

        // register all weka object editors (TODO: in headless mode actually not
        // necessary)
        GenericObjectEditor.registerEditors();

        // load all weka classes from the weka library
        try {
            Properties props = new Properties();
            props.load(new Utils().getClass().getClassLoader()
                    .getResources("weka/gui/GenericObjectEditor.props")
                    .nextElement().openStream());
            PluginManager.addFromProperties(props);

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private GenericObjectEditor m_ClassifierEditor;

    /**
     * @param model
     */
    public DialogComponentWekaClassifier(SettingsModelWekaClassifier model,
            String label) {
        super(model);

        // Add Weka panel for selecting the classifier and its options
        m_ClassifierEditor = new GenericObjectEditor();
        PropertyPanel m_CEPanel = new PropertyPanel(m_ClassifierEditor);
        m_ClassifierEditor.setClassType(Classifier.class);

        JPanel panel = new JPanel(new FlowLayout());
        panel.add(new JLabel(label));
        panel.add(m_CEPanel);

        // add classifier editor panel
        getComponentPanel().add(panel);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void updateComponent() {
        m_ClassifierEditor.setValue(((SettingsModelWekaClassifier)getModel())
                .getClassifier());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettingsBeforeSave() throws InvalidSettingsException {
        // update model
        ((SettingsModelWekaClassifier)getModel())
                .setClassifier((AbstractClassifier)m_ClassifierEditor
                        .getValue());

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void checkConfigurabilityBeforeLoad(PortObjectSpec[] specs)
            throws NotConfigurableException {
        // TODO Auto-generated method stub

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void setEnabledComponents(boolean enabled) {
        m_ClassifierEditor.setEnabled(enabled);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setToolTipText(String text) {
        // nothing to do
    }

}
