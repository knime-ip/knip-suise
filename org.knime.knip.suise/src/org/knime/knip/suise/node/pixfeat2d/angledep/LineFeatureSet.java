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

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

import org.knime.knip.suise.data.feat.Dart;
import org.knime.knip.suise.data.feat.Dart2D;
import org.knime.knip.suise.data.feat.FeatImgGenerator;
import org.knime.knip.suise.data.feat.PixFeatureSet;

/**
 * 
 * @param <T> image type
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class LineFeatureSet<T extends RealType<T>, F extends RealType<F>>
        implements PixFeatureSet<T>, FeatImgGenerator<F> {

    private int m_radius;

    private String[] m_names;

    private RandomAccess<T> m_srcRA;

    private Img<T> m_srcImg;

    private long[] m_pos;

    private long[] m_tmpPos = new long[2];

    private long m_tmp;

    private double m_angle;

    /**
     * @param radius
     * @param target
     */
    public LineFeatureSet(final int radius) {
        super();

        m_radius = radius;
        init();

    }

    private void init() {
        m_names = new String[m_radius * 2 + 1];
        for (int i = 0; i < m_radius * 2 + 1; i++) {
            m_names[i] = "LineFeature[pos=" + i + "]";
        }
    }

    public LineFeatureSet() {
        // no-op
    }

    // /**
    // * {@inheritDoc}
    // */
    // @Override
    // public void assignFeatureValues(LocalizablePlaneCursor<FloatType> cursor,
    // int featureID, double angle) {
    //
    // int angleID = (int) (Math.round(angle / (2 * Math.PI) * 2 * m_numAng));
    // T minVal = m_img.createType();
    // minVal.setReal(minVal.getMinValue());
    // LocalizableByDimCursor<T> srcCursor = m_img
    // .createLocalizableByDimCursor(new OutOfBoundsStrategyValueFactory<T>(
    // minVal));
    // srcCursor.fwd();
    // int x = (int) Math.round(Math.sin(angle) * (featureID - m_radius));
    // int y = (int) Math.round(Math.cos(angle) * (featureID - m_radius));
    //
    // // // garanties that the manhatten-distance between the two IntPoints is
    // // // twice
    // // // the radius (instead of the euclidean distance)
    // // if (Math.abs(x) == Math.abs(y) && Math.abs(x) != m_radius) {
    // // x = m_radius * x / Math.abs(x);
    // // y = m_radius * y / Math.abs(y);
    // // }
    // // x = (Math.abs(x) > Math.abs(y) && Math.abs(x) != m_radius) ? m_radius
    // // * x / Math.abs(x) : x;
    // // y = (Math.abs(y) > Math.abs(x) && Math.abs(y) != m_radius) ? m_radius
    // // * y / Math.abs(y) : y;
    // int[] offset = new int[] { -y, x };
    // int[] pos = new int[2];
    //
    // while (cursor.hasNext()) {
    // cursor.fwd();
    // pos[0] = cursor.getPosition(0);
    // pos[1] = cursor.getPosition(1);
    // srcCursor.setPosition(pos[0] + offset[0], 0);
    // srcCursor.setPosition(pos[1] + offset[1], 1);
    // cursor.getType().set((srcCursor.getType().getRealFloat()));
    // }
    // srcCursor.close();
    //
    // }

    /**
     * {@inheritDoc}
     */
    @Override
    public double value(int id) {
        long x = m_pos[0];
        long y = m_pos[1];
        m_tmpPos[0] = (long)Math.round(Math.sin(m_angle) * (id - m_radius));
        m_tmpPos[1] = (long)Math.round(Math.cos(m_angle) * (id - m_radius));

        // // garanties that the manhatten-distance between the two
        // IntPoints is
        // // twice
        // // the radius (instead of the euclidean distance)
        // if (Math.abs(x) == Math.abs(y) && Math.abs(x) != m_radius) {
        // x = m_radius * x / Math.abs(x);
        // y = m_radius * y / Math.abs(y);
        // }
        // x = (Math.abs(x) > Math.abs(y) && Math.abs(x) != m_radius) ?
        // m_radius
        // * x / Math.abs(x) : x;
        // y = (Math.abs(y) > Math.abs(x) && Math.abs(y) != m_radius) ?
        // m_radius
        // * y / Math.abs(y) : y;

        // orthogonal vector
        m_tmp = m_tmpPos[0];
        m_tmpPos[0] = -m_tmpPos[1];
        m_tmpPos[1] = m_tmp;

        m_tmpPos[0] += x;
        m_tmpPos[1] += y;

        m_srcRA.setPosition(m_tmpPos[0], 0);
        m_srcRA.setPosition(m_tmpPos[1], 1);
        return m_srcRA.get().getRealDouble();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numFeatures() {
        return m_names.length;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String featureSetId() {
        return "Line Feature Factory";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void updateImg(Img<T> img) {
        T minVal = img.firstElement().createVariable();
        minVal.setReal(minVal.getMinValue());
        m_srcRA = Views.extendValue(img, minVal).randomAccess();
        m_srcImg = img;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void updateDart(Dart dart) {
        if (dart instanceof Dart2D) {
            m_angle = ((Dart2D)dart).angle();
            m_pos = dart.offset();
        } else {
            throw new IllegalArgumentException(
                    "3D not supported yet in the line feature set.");
        }

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void generateFeatImg(RandomAccessibleInterval<F> featImg,
            int featID, Dart dart) {

        if (dart instanceof Dart2D) {
            m_angle = ((Dart2D)dart).angle();
        } else {
            throw new IllegalArgumentException(
                    "3D not supported yet in the line feature set.");
        }

        m_tmpPos[0] = (long)Math.round(Math.sin(m_angle) * (featID - m_radius));
        m_tmpPos[1] = (long)Math.round(Math.cos(m_angle) * (featID - m_radius));

        // orthogonal vector
        m_tmp = m_tmpPos[0];
        m_tmpPos[0] = -m_tmpPos[1];
        m_tmpPos[1] = m_tmp;

        IterableInterval<T> translated =
                Views.iterable(Views.offset(m_srcImg, m_tmpPos));

        Cursor<T> srcC = translated.localizingCursor();
        RandomAccess<F> resRA = Views.extendBorder(featImg).randomAccess();
        while (srcC.hasNext()) {
            srcC.fwd();
            resRA.setPosition(srcC);
            resRA.get().setReal(srcC.get().getRealDouble());
        }

    }

    // /**
    // * {@inheritDoc}
    // */
    // @Override
    // public void writeExternal(ObjectOutput out) throws IOException {
    // super.writeExternal(out);
    // out.writeInt(m_radius);
    // }
    //
    // /**
    // * {@inheritDoc}
    // */
    // @Override
    // public void readExternal(ObjectInput in) throws IOException,
    // ClassNotFoundException {
    // super.readExternal(in);
    // m_radius = in.readInt();
    // init();
    //
    // }

}
