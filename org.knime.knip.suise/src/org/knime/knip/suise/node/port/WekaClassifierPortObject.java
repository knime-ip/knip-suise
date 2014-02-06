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

import java.awt.GridLayout;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.zip.ZipEntry;

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import org.apache.mahout.math.Arrays;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.util.NonClosableInputStream;
import org.knime.core.data.util.NonClosableOutputStream;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.ModelContentRO;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.PortObjectSpec;
import org.knime.core.node.port.PortObjectZipInputStream;
import org.knime.core.node.port.PortObjectZipOutputStream;
import org.knime.core.node.port.PortType;

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.Classifier;
import weka.classifiers.trees.j48.Distribution;
import weka.core.Instance;
import weka.core.Instances;

/**
 * 
 * @author hornm, University of Konstanz
 */
public class WekaClassifierPortObject implements PortObject {

    /**
     * The Port Type.
     */
    public static final PortType TYPE = new PortType(
            WekaClassifierPortObject.class);

    private static final NodeLogger LOGGER = NodeLogger
            .getLogger(WekaClassifierPortObject.class);

    /**
     * @return Serializer for the {@link WekaClassifierModelPortObject}
     */
    public static PortObjectSerializer<WekaClassifierPortObject> getPortObjectSerializer() {
        return new PortObjectSerializer<WekaClassifierPortObject>() {

            /** {@inheritDoc} */
            @Override
            public void savePortObject(
                    final WekaClassifierPortObject portObject,
                    final PortObjectZipOutputStream out,
                    final ExecutionMonitor exec) throws IOException,
                    CanceledExecutionException {
                portObject.save(out);

            }

            /** {@inheritDoc} */
            @Override
            public WekaClassifierPortObject loadPortObject(
                    final PortObjectZipInputStream in,
                    final PortObjectSpec spec, final ExecutionMonitor exec)
                    throws IOException, CanceledExecutionException {
                return load(in, (WekaClassifierPortObjectSpec)spec);
            }
        };
    }

    private void save(final PortObjectZipOutputStream out) {
        // save weka classifier
        ObjectOutputStream oo = null;
        try {
            out.putNextEntry(new ZipEntry("classifier.objectout"));
            oo = new ObjectOutputStream(new NonClosableOutputStream.Zip(out));
            oo.writeObject(m_classifier);
        } catch (IOException ioe) {
            LOGGER.error("Internal error: Could not save settings", ioe);
        } finally {
            if (oo != null) {
                try {
                    oo.close();
                } catch (Exception e) {
                    LOGGER.debug("Could not close stream", e);
                }
            }
        }

        // save training instances
        try {
            out.putNextEntry(new ZipEntry("training.objectout"));
            oo = new ObjectOutputStream(new NonClosableOutputStream.Zip(out));

            // only write the instance metadata, not the instances itself
            oo.writeObject(new Instances(m_trainingInstances, 0));
        } catch (IOException ioe) {
            LOGGER.error("Internal error: Could not save settings", ioe);
        } finally {
            if (oo != null) {
                try {
                    oo.close();
                } catch (Exception e) {
                    LOGGER.debug("Could not close stream", e);
                }
            }
        }

        // // save meta information
        // ModelContent model = new ModelContent(MODEL_INFO);
        // Config mapperconf = model.addConfig(MAPPER_KEY);
        // m_mapper.save(mapperconf);
        // try {
        // out.putNextEntry(new ZipEntry("mapper.xmlout"));
        // model.saveToXML(out);
        // } catch (IOException ioe) {
        // LOGGER.error("Internal error: Could not save settings", ioe);
        // }
    }

    private static WekaClassifierPortObject load(
            final PortObjectZipInputStream in,
            final WekaClassifierPortObjectSpec spec) {
        ObjectInputStream oi = null;
        Classifier classifier = null;
        Instances trainInstances = null;
        ModelContentRO model = null;

        try {
            // load classifier
            ZipEntry zentry = in.getNextEntry();
            assert zentry.getName().equals("classifier.objectout");
            oi = new ObjectInputStream(new NonClosableInputStream.Zip(in));
            classifier = (Classifier)oi.readObject();
        } catch (IOException ioe) {
            LOGGER.error("Internal error: Could not load settings", ioe);
        } catch (ClassNotFoundException cnf) {
            LOGGER.error("Internal error: Could not load settings", cnf);
        } finally {
            if (oi != null) {
                try {
                    oi.close();
                } catch (Exception e) {
                    LOGGER.debug("Could not close stream", e);
                }
            }
        }

        try {
            // load training instances
            ZipEntry zentry = in.getNextEntry();
            assert zentry.getName().equals("training.objectout");
            oi = new ObjectInputStream(new NonClosableInputStream.Zip(in));
            trainInstances = (Instances)oi.readObject();
        } catch (IOException ioe) {
            LOGGER.error("Internal error: Could not load settings", ioe);
        } catch (ClassNotFoundException cnf) {
            LOGGER.error("Internal error: Could not load settings", cnf);
        } finally {
            if (oi != null) {
                try {
                    oi.close();
                } catch (Exception e) {
                    LOGGER.debug("Could not close stream", e);
                }
            }
        }

        // try {
        // // load meta info
        // ZipEntry zentry = in.getNextEntry();
        // assert zentry.getName().equals("mapper.xmlout");
        // model = ModelContent.loadFromXML(in);
        // } catch (IOException ioe) {
        // LOGGER.error("Internal error: Could not load settings", ioe);
        // }
        assert (classifier != null);
        assert (trainInstances != null);
        assert (model != null);

        // DataCellStringMapper mapper = null;
        // try {
        // mapper = DataCellStringMapper.load(model.getConfig(MAPPER_KEY));
        // } catch (InvalidSettingsException ise) {
        // LOGGER.error("Internal error: Could not load settings", ise);
        // }
        return new WekaClassifierPortObject(classifier, trainInstances, spec);
    }

    private Classifier m_classifier;

    private Instances m_trainingInstances;

    // private DataCellStringMapper m_mapper;

    private WekaClassifierPortObjectSpec m_modelspec;

    /**
     * The WekaClassifierPortObject holds information about the used classifier,
     * training instances, columns and class column.
     * 
     * @param classifier Classifier from weka.
     * @param traininginstances training instances used.
     * @param spec the {@link WekaClassifierModelPortObjectSpec}.
     */
    public WekaClassifierPortObject(final Classifier classifier,
            final Instances traininginstances,
            final WekaClassifierPortObjectSpec spec) {
        m_classifier = classifier;
        // m_trainingInstances = new Instances(traininginstances, 0);
        m_trainingInstances = traininginstances;
        m_modelspec = spec;
    }

    /**
     * @return the classifier
     */
    public Classifier getClassifier() {
        return m_classifier;
    }

    /**
     * @return the training {@link DataTableSpec}.
     */
    public WekaClassifierPortObjectSpec getSpec() {
        return m_modelspec;
    }

    /** {@inheritDoc} */
    @Override
    public String getSummary() {
        return m_classifier.getClass().getSimpleName();
    }

    // /**
    // * @return the mapper, mapping DataCells to Strings.
    // */
    // public DataCellStringMapper getMapper() {
    // return m_mapper;
    // }

    /**
     * @return the training Instances
     */
    public Instances getTrainingInstances() {
        return m_trainingInstances;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public JComponent[] getViews() {

        JPanel panel = new JPanel(new GridLayout(1, 2));

        // classifier info
        JPanel classifierInfoPanel = new JPanel();
        String text = m_classifier.getClass().getSimpleName() + "\n";
        text += m_classifier.toString();
        classifierInfoPanel.add(new JTextArea(text));
        classifierInfoPanel.setBorder(BorderFactory
                .createTitledBorder("Weka Classifier Details"));
        panel.add(new JScrollPane(classifierInfoPanel));

        // training instances info
        if (m_trainingInstances.numInstances() > 0) {
            JPanel dataInfoPanel = new JPanel();
            text = "Training instances\n\n";
            Distribution distr =
                    new Distribution(1, m_trainingInstances.numClasses());
            for (Instance i : m_trainingInstances) {
                try {
                    distr.add(0, i);
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }

            text +=
                    "Class-distribution: " + Arrays.toString(distr.matrix()[0])
                            + " (labels: "
                            + Arrays.toString(m_modelspec.getClassLabels())
                            + ")" + "\n\n";
            text += m_trainingInstances.toSummaryString();

            dataInfoPanel.add(new JTextArea(text));
            dataInfoPanel.setBorder(BorderFactory
                    .createTitledBorder("Training Instances Details"));
            panel.add(new JScrollPane(dataInfoPanel));
        }

        panel.setName("Weka Outport");
        return new JComponent[]{panel};
    }
}
