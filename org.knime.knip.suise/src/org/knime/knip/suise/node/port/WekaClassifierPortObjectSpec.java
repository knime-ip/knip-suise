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
package org.knime.knip.suise.node.port;

import java.io.IOException;
import java.util.zip.ZipEntry;

import javax.swing.JComponent;

import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.ModelContent;
import org.knime.core.node.ModelContentRO;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.port.PortObjectSpec;
import org.knime.core.node.port.PortObjectSpecZipInputStream;
import org.knime.core.node.port.PortObjectSpecZipOutputStream;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class WekaClassifierPortObjectSpec implements PortObjectSpec {

    private static final NodeLogger LOGGER = NodeLogger
            .getLogger(WekaClassifierPortObjectSpec.class);

    /*
     * Model info identifier.
     */
    private static final String MODEL_INFO = "model_info";

    /*
     * Key to store the class labels.
     */
    private static final String CLASSLABELS_KEY = "class_labels";

    /**
     * @return Serializer for the {@link WekaClassifierPortObjectSpec}.
     */
    public static PortObjectSpecSerializer<WekaClassifierPortObjectSpec> getPortObjectSpecSerializer() {
        return new PortObjectSpecSerializer<WekaClassifierPortObjectSpec>() {
            /** {@inheritDoc} */
            @Override
            public void savePortObjectSpec(
                    final WekaClassifierPortObjectSpec portObject,
                    final PortObjectSpecZipOutputStream out) throws IOException {
                portObject.save(out);

            }

            /** {@inheritDoc} */
            @Override
            public WekaClassifierPortObjectSpec loadPortObjectSpec(
                    final PortObjectSpecZipInputStream in) throws IOException {
                return load(in);
            }
        };
    }

    private void save(final PortObjectSpecZipOutputStream out) {
        ModelContent model = new ModelContent(MODEL_INFO);
        model.addStringArray(CLASSLABELS_KEY, m_classLabels);
        try {
            out.putNextEntry(new ZipEntry("classlabels.xmlout"));
            model.saveToXML(out);
        } catch (IOException ioe) {
            LOGGER.error("Internal error: Could not save settings", ioe);
        }
    }

    private static WekaClassifierPortObjectSpec load(
            final PortObjectSpecZipInputStream in) {
        ModelContentRO model = null;
        try {
            ZipEntry zentry = in.getNextEntry();
            assert zentry.getName().equals("classlabels.xmlout");
            model = ModelContent.loadFromXML(in);
        } catch (IOException ioe) {
            LOGGER.error("Internal error: Could not load settings", ioe);
        }
        String[] classLabels = null;
        try {
            classLabels = model.getStringArray(CLASSLABELS_KEY);
        } catch (InvalidSettingsException ise) {
            LOGGER.error("Internal error: Could not load settings", ise);
        }
        return new WekaClassifierPortObjectSpec(classLabels);
    }

    private final String[] m_classLabels;

    /**
     * The {@link WekaClassifierPortObjectSpec} holds the columns of the
     * training data and the class column.
     * 
     * @param classLabels the class labels used to train the classifier
     * 
     */
    public WekaClassifierPortObjectSpec(String[] classLabels) {
        m_classLabels = classLabels;
    }

    /**
     * @return the class labels used to train the classifier
     */
    public String[] getClassLabels() {
        return m_classLabels;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public JComponent[] getViews() {
        return new JComponent[]{};
    }
}
