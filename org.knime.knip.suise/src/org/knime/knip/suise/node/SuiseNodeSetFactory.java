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
package org.knime.knip.suise.node;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSetFactory;
import org.knime.core.node.config.ConfigRO;
import org.knime.knip.suise.node.boundarymodel.BoundaryModelNodeFactory;
import org.knime.knip.suise.node.levelcca.LevelCCANodeFactory;
/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import org.knime.knip.suise.node.pixclassapply.PixClassApplyNodeFactory;
import org.knime.knip.suise.node.pixclassmodel.PixClassModelNodeFactory;
import org.knime.knip.suise.node.pixfeat2d.angledep.AngleDepPixFeat2DNodeFactory;
import org.knime.knip.suise.node.pixfeat2d.fiji.FijiTrainSegFeaturesNodeFactory;

import weka.core.Environment;
import weka.gui.GenericObjectEditor;

/**
 * 
 * @author Martin Horn, University of Konstanz
 */
public class SuiseNodeSetFactory implements NodeSetFactory {

    static {
        // deactivate the package management of weka
        Environment env = Environment.getSystemWide();
        env.addVariableSystemWide("weka.core.loadPackages", "false");
        env.addVariableSystemWide("weka.packageManager.offline", "true");

        // register all weka object editors (TODO: in headless mode actually not
        // necessary)
        GenericObjectEditor.registerEditors();

        // make the classes available within weka
        // try {
        // // add the classifieres etc. from the weka-classes.props-file
        // PluginManager.addFromProperties(wekaProperties);
        // } catch (Exception e) {
        // throw new RuntimeException("Could not initialize weka nodes!", e);
        // }
    }

    private final Map<String, String> m_nodeFactories =
            new HashMap<String, String>();

    /**
     * {@inheritDoc}
     */
    @Override
    public ConfigRO getAdditionalSettings(final String id) {
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getAfterID(final String id) {
        return "";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getCategoryPath(final String id) {
        return m_nodeFactories.get(id);
    }

    /**
     * {@inheritDoc}
     */
    @SuppressWarnings("unchecked")
    @Override
    public Class<? extends NodeFactory<? extends NodeModel>> getNodeFactory(
            final String id) {
        try {
            return (Class<? extends NodeFactory<? extends NodeModel>>)Class
                    .forName(id);
        } catch (final ClassNotFoundException e) {
        }
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Collection<String> getNodeFactoryIds() {
        m_nodeFactories.put(
                FijiTrainSegFeaturesNodeFactory.class.getCanonicalName(),
                "/community/knip/suise");
        m_nodeFactories.put(PixClassModelNodeFactory.class.getCanonicalName(),
                "/community/knip/suise");
        m_nodeFactories.put(PixClassApplyNodeFactory.class.getCanonicalName(),
                "/community/knip/suise");
        m_nodeFactories.put(BoundaryModelNodeFactory.class.getCanonicalName(),
                "/community/knip/suise");
        m_nodeFactories.put(
                AngleDepPixFeat2DNodeFactory.class.getCanonicalName(),
                "/community/knip/suise");
        // m_nodeFactories.put(
        // LabelingCompareNodeFactory.class.getCanonicalName(),
        // "/community/knip/suise");
        m_nodeFactories.put(LevelCCANodeFactory.class.getCanonicalName(),
                "/community/knip/suise");

        return m_nodeFactories.keySet();
    }
}
