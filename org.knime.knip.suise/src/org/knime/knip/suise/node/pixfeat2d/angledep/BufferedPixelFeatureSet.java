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

import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

import org.knime.knip.suise.data.feat.Dart;
import org.knime.knip.suise.data.feat.PixFeatureSet;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class BufferedPixelFeatureSet<T extends RealType<T>> implements
        PixFeatureSet {

    private RandomAccess<T> m_featureImageAccess;

    private int m_currentAngleID;

    /* fields to be made persistent to restore the factory */
    private String[] m_featureNames;

    private int m_xDim;

    private int m_yDim;

    private int m_featureDimension;

    private int m_angularDimension;

    private long[] m_position;

    private Img m_img;

    /**
     * For serialization purposes. Must not be called by the user.
     */
    public BufferedPixelFeatureSet() {
        // no-op
    }

    public BufferedPixelFeatureSet(String[] featureNames, int xDim, int yDim,
            int featureDimension) {
        this(featureNames, xDim, yDim, featureDimension, -1);

    }

    public BufferedPixelFeatureSet(String[] featureNames, int xDim, int yDim,
            int featureDimension, int angularDimension) {
        m_featureNames = featureNames;
        m_xDim = xDim;
        m_yDim = yDim;
        m_featureDimension = featureDimension;
        m_angularDimension = angularDimension;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double value(int id) {
        if (m_angularDimension != -1) {
            m_position[m_angularDimension] = m_currentAngleID;
        }
        m_position[m_featureDimension] = id;
        m_featureImageAccess.setPosition(m_position);
        return m_featureImageAccess.get().getRealDouble();

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numFeatures() {
        return m_featureNames.length;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String featureSetId() {
        return "Buffered Pixel Feature Factory";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void updateImg(Img img) {
        if (img.numDimensions() < m_xDim || img.numDimensions() < m_yDim
                || img.numDimensions() < m_featureDimension
                || img.numDimensions() < m_angularDimension) {
            throw new IllegalArgumentException(
                    "Wrong number of dimension in image.");
        }
        if (img.dimension(m_featureDimension) != m_featureNames.length) {
            throw new IllegalArgumentException(
                    "The feature dimension of the image doesn't fit the number of features.");
        }
        // if (m_angularDimension != -1) {
        // m_currentAngleID = (int) Math.round((m_angle / (2 * Math.PI))
        // * (double) m_img.dimension(m_angularDimension))
        // % (int) m_img.dimension(m_angularDimension);
        // }

        m_position = new long[img.numDimensions()];

        m_featureImageAccess = Views.extendBorder(img).randomAccess();

        m_img = img;

    }

    /**
     * {@inheritDoc}
     */
    public void updateDart(Dart dart) {
        if (m_angularDimension != -1) {
            // m_currentAngleID = (int) Math.round((m_angle / (2 * Math.PI))
            // * (double) m_img.dimension(m_angularDimension))
            // % (int) m_img.dimension(m_angularDimension);
            m_currentAngleID = dart.directionIndex();
        }

        m_position[m_xDim] = dart.offset()[0];
        m_position[m_yDim] = dart.offset()[1];

    }

    public int numAngles() {
        if (m_img == null) {
            throw new IllegalStateException("Image not set!");
        } else {
            return (int)m_img.dimension(m_angularDimension);
        }
    }

    public BufferedPixelFeatureSet<T> copy() {
        BufferedPixelFeatureSet<T> res =
                new BufferedPixelFeatureSet<T>(m_featureNames, m_xDim, m_yDim,
                        m_featureDimension, m_angularDimension);
        res.updateImg(m_img);
        return res;
    }
    // /**
    // * {@inheritDoc}
    // */
    // @Override
    // public void writeExternal(ObjectOutput out) throws IOException {
    // super.writeExternal(out);
    // out.writeInt(getNumTotalFeatures());
    // for (int i = 0; i < getNumTotalFeatures(); i++) {
    // out.writeUTF(m_featureNames[i]);
    // }
    // out.writeInt(m_xDim);
    // out.writeInt(m_yDim);
    // out.writeInt(m_featureDimension);
    // out.writeInt(m_angularDimension);
    //
    // }

    // /**
    // * {@inheritDoc}
    // */
    // @Override
    // public void readExternal(ObjectInput in) throws IOException,
    // ClassNotFoundException {
    // super.readExternal(in);
    // m_featureNames = new String[in.readInt()];
    // for (int i = 0; i < m_featureNames.length; i++) {
    // m_featureNames[i] = in.readUTF();
    // }
    // m_xDim = in.readInt();
    // m_yDim = in.readInt();
    // m_featureDimension = in.readInt();
    // m_angularDimension = in.readInt();
    // }

}
