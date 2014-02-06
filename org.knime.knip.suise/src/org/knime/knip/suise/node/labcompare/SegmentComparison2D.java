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
package org.knime.knip.suise.node.labcompare;

import java.util.List;
import java.util.Vector;

import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.ops.operation.randomaccessibleinterval.unary.DistanceMap;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

import org.knime.core.node.NodeLogger;
import org.knime.knip.core.features.seg.ExtractOutlineImg;
import org.knime.knip.core.ops.img.IterableIntervalNormalize;

/**
 * Compares two image segments (represented as binary image) usually from
 * different images including different measures (see article Lin, G.; K., A.;
 * Olson, K.; Guzowski, J. F.; Barnes, C. A. & Roysam, B. A hybrid 3d watershed
 * algorithm incorporating gradient cues and object models for automatic
 * segmentation of nuclei in confocal image stacks Cytometry Part A, 2003, 56A,
 * 23-36)
 * 
 * - the hausdorff metric
 * 
 * - number of false positive and negative pixels
 * 
 * - number of overlapping pixels
 * 
 * - NSD
 * 
 * 
 * 
 * 
 */

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class SegmentComparison2D {

    private int m_fp; // false positive

    private int m_fn; // false negative

    private int m_total;

    private double m_hausdorffdist;

    private double m_nsd;

    // private double m_nsd; // normalised sum of distances

    private Img<BitType> m_seg1;

    private Img<BitType> m_seg2;

    private Img<ByteType> m_overlayImg;

    private static final byte SEG1_VAL = -100;

    private static final byte SEG2_VAL = -75;

    private static final byte SEG1_SEG2_VAL = -77;

    private static final byte BG = 0;

    /**
     * 
     * 
     * @param seg1 - segment from the reference image
     * @param seg1Pos its position in the source image
     * @param seg2 - segment from the segmented image
     * @param seg2Pos its position in the source image
     */
    public SegmentComparison2D(final Segment seg1, final Segment seg2) {

        if (seg1.mask.numDimensions() > 2 || seg2.mask.numDimensions() > 2) {
            NodeLogger
                    .getLogger(LabelingComparison2D.class)
                    .warn("Segment comparision only supported for 2D-segments. Dimension 0 and 1 are chosen.");
        }

        long x = Math.min(seg1.offset[0], seg2.offset[0]);
        long y = Math.min(seg1.offset[1], seg2.offset[1]);

        long width =
                Math.max(seg1.offset[0] + seg1.mask.dimension(0),
                        seg2.offset[0] + seg2.mask.dimension(0))
                        - Math.min(seg1.offset[0], seg2.offset[0]);
        long height =
                Math.max(seg1.offset[1] + seg1.mask.dimension(1),
                        seg2.offset[1] + seg2.mask.dimension(1))
                        - Math.min(seg1.offset[1], seg2.offset[1]);

        ImgFactory<BitType> imgFac = new ArrayImgFactory<BitType>();
        m_seg1 = imgFac.create(new long[]{width, height}, new BitType());
        m_seg2 = imgFac.create(new long[]{width, height}, new BitType());

        Cursor<BitType> newSeg1Cur = m_seg1.localizingCursor();
        RandomAccess<BitType> seg1RA =
                Views.extendValue(seg1.mask, new BitType(false)).randomAccess();

        while (newSeg1Cur.hasNext()) {
            newSeg1Cur.fwd();
            seg1RA.setPosition(newSeg1Cur.getIntPosition(0) - seg1.offset[0]
                    + x, 0);
            seg1RA.setPosition(newSeg1Cur.getIntPosition(1) - seg1.offset[1]
                    + y, 1);
            newSeg1Cur.get().set(seg1RA.get());
        }

        Cursor<BitType> newSeg2Cur = m_seg2.localizingCursor();
        RandomAccess<BitType> seg2RA =
                Views.extendValue(seg2.mask, new BitType(false)).randomAccess();

        while (newSeg2Cur.hasNext()) {
            newSeg2Cur.fwd();
            seg2RA.setPosition(newSeg2Cur.getIntPosition(0) - seg2.offset[0]
                    + x, 0);
            seg2RA.setPosition(newSeg2Cur.getIntPosition(1) - seg2.offset[1]
                    + y, 1);
            newSeg2Cur.get().set(seg2RA.get());
        }

        m_overlayImg =
                new ArrayImgFactory<ByteType>().create(
                        new long[]{width, height}, new ByteType());

        RandomAccess<ByteType> overlayRA = m_overlayImg.randomAccess();

        // m_ref.showInFrame();
        // m_seg.showInFrame();

        m_fn = m_fp = m_total = 0;

        newSeg1Cur.reset();
        newSeg2Cur.reset();
        BitType seg2Val = newSeg2Cur.get();
        BitType seg1Val = newSeg1Cur.get();

        while (newSeg1Cur.hasNext()) {
            newSeg1Cur.fwd();
            newSeg2Cur.fwd();

            if (seg1Val.get() && !seg2Val.get()) {
                m_fn++;
                overlayRA.setPosition(newSeg1Cur);
                overlayRA.get().set(SEG1_VAL);
            } else if (!seg1Val.get() && seg2Val.get()) {
                m_fp++;
                overlayRA.setPosition(newSeg1Cur);
                overlayRA.get().set(SEG2_VAL);
            }
            if (seg1Val.get() || seg2Val.get()) {
                m_total++;

            }

            if (seg1Val.get() && seg2Val.get()) {
                overlayRA.setPosition(newSeg1Cur);
                overlayRA.get().set(SEG1_SEG2_VAL);
            }

        }

        calcNSD();
        calcHausdorffDist();

    }

    /**
     * @return number of false positive pixels
     */

    public int getFPPixels() {
        return m_fp;
    }

    /**
     * @return number of false negative pixels
     */

    public int getFNPixels() {
        return m_fn;
    }

    public int getTotalPixels() {
        return m_total;
    }

    /**
     * @return the hausdorff distance
     */

    public double getHausdorffDistance() {
        return m_hausdorffdist;
    }

    /**
     * 
     * @return the normalized sum of distances
     */
    public double getNSD() {
        return m_nsd;
    }

    private void calcHausdorffDist() {

        // brute force

        // get border pixels
        List<int[]> borderSeg1 = new Vector<int[]>();
        List<int[]> borderSeg2 = new Vector<int[]>();

        Cursor<BitType> seg1Cur = m_seg1.localizingCursor();
        Cursor<BitType> seg2Cur = m_seg2.localizingCursor();

        RandomAccess<BitType> seg1NRA =
                Views.extendValue(m_seg1, new BitType(false)).randomAccess();
        RandomAccess<BitType> seg2NRA =
                Views.extendValue(m_seg2, new BitType(false)).randomAccess();

        // RandomAccess<ByteType> debugRA = m_overlayImg.randomAccess();
        // Cursor<ByteType> debugCur = m_overlayImg.cursor();

        while (seg1Cur.hasNext()) {
            seg1Cur.fwd();
            seg2Cur.fwd();
            // debugCur.fwd();

            if (seg1Cur.get().get() && numNeighborPixels(seg1NRA, seg1Cur) < 8) {
                int[] pos = new int[seg1Cur.numDimensions()];
                seg1Cur.localize(pos);
                borderSeg1.add(pos);
                // debugCur.get().set((byte)0);
            }
            if (seg2Cur.get().get() && numNeighborPixels(seg2NRA, seg2Cur) < 8) {
                int[] pos = new int[seg2Cur.numDimensions()];
                seg2Cur.localize(pos);
                borderSeg2.add(pos);
                // debugCur.get().set((byte)0);
            }
        }

        double[][] distMatrix = distMatrix(borderSeg1, borderSeg2);
        double maxDist =
                Math.sqrt(Math.pow(m_seg1.dimension(0), 2)
                        + Math.pow(m_seg1.dimension(1), 2));

        Point[] maxP1 = new Point[]{new Point(2), new Point(2)};
        Point[] maxP2 = new Point[]{new Point(2), new Point(2)};
        Point tmp = new Point(2);

        // h(ref,seq)
        double hrs = 0;
        for (int i = 0; i < borderSeg1.size(); i++) {
            double min = maxDist;
            for (int j = 0; j < borderSeg2.size(); j++) {

                if (min < distMatrix[i][j]) {
                    tmp = new Point(borderSeg2.get(j));
                }

                min = Math.min(min, distMatrix[i][j]);
            }

            if (min > hrs) {
                maxP1[0] = tmp;
                maxP1[1] = new Point(borderSeg1.get(i));
            }

            hrs = Math.max(hrs, min);
        }

        // h(seq,ref)
        double hsr = 0;
        for (int i = 0; i < borderSeg2.size(); i++) {
            double min = maxDist;
            for (int j = 0; j < borderSeg1.size(); j++) {

                if (min < distMatrix[j][i]) {
                    tmp = new Point(borderSeg1.get(j));
                }

                min = Math.min(min, distMatrix[j][i]);
            }
            if (min > hsr) {
                maxP2[0] = tmp;
                maxP2[1] = new Point(borderSeg2.get(i));
            }

            hsr = Math.max(hsr, min);
        }

        // debugRA.setPosition(maxP1[0]);
        // debugRA.get().set((byte)128);
        // debugRA.setPosition(maxP1[1]);
        // debugRA.get().set((byte)128);
        // debugRA.setPosition(maxP2[0]);
        // debugRA.get().set((byte)128);
        // debugRA.setPosition(maxP2[1]);
        // debugRA.get().set((byte)128);

        m_hausdorffdist = Math.max(hrs, hsr);
        // m_debugImage.setTitle("overlap: "
        // + (1.0 - (double) (m_fn + m_fp) / (double) m_total));
        // m_debugImage.setTitle("hd: " + m_hausdorffdist);

    }

    private void calcNSD() {
        // create outline image from the reference segment within the overlay
        // image
        Img<BitType> outline =
                new ArrayImgFactory<BitType>().create(m_overlayImg,
                        new BitType());
        Cursor<BitType> outlineCur = outline.cursor();
        Cursor<ByteType> overlayCur = m_overlayImg.cursor();
        while (outlineCur.hasNext()) {
            outlineCur.fwd();
            overlayCur.fwd();
            if (overlayCur.get().get() == SEG1_VAL
                    || overlayCur.get().get() == SEG1_SEG2_VAL) {
                outlineCur.get().set(true);
            }
        }
        Img<BitType> tmp =
                new ArrayImgFactory<BitType>().create(m_overlayImg,
                        new BitType());
        new ExtractOutlineImg(false).compute(outline, tmp);
        outline = tmp;
        // invert outline
        for (BitType type : outline) {
            type.set(!type.get());
        }
        Img<FloatType> distTrans =
                new ArrayImgFactory<FloatType>().create(outline,
                        new FloatType());
        new DistanceMap<BitType>().compute(outline, distTrans);
        new IterableIntervalNormalize<FloatType>(0, new FloatType(), null, true)
                .compute(distTrans, distTrans);

        // sum distance over the union of all pixels and non-overlapping pixels
        overlayCur.reset();
        Cursor<FloatType> distCur = distTrans.cursor();
        double union = 0;
        double no_overlap = 0;
        while (overlayCur.hasNext()) {
            overlayCur.fwd();
            distCur.fwd();
            if (overlayCur.get().get() == SEG1_VAL
                    || overlayCur.get().get() == SEG2_VAL) {
                no_overlap += distCur.get().get();
            }
            if (overlayCur.get().get() != BG) {
                union += distCur.get().get();
            }
        }
        m_nsd = no_overlap / union;

    }

    private int numNeighborPixels(RandomAccess<BitType> nRA,
            Cursor<BitType> orgCur) {
        // TODO: optimization with LocalNeighborhoodCursor
        BitType currentVal = orgCur.get();
        // center position
        int[] c = new int[orgCur.numDimensions()];
        orgCur.localize(c);
        nRA.setPosition(orgCur);
        int num = 0;
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                nRA.setPosition(c[0] + i, 0);
                nRA.setPosition(c[1] + j, 1);
                if (nRA.get().compareTo(currentVal) == 0) {
                    num++;
                }
            }
        }

        return num;

    }

    private double[][] distMatrix(final List<int[]> v1, final List<int[]> v2) {

        double[][] m = new double[v1.size()][v2.size()];
        for (int i = 0; i < v1.size(); i++) {
            for (int j = 0; j < v2.size(); j++) {
                m[i][j] =
                        Math.sqrt(Math.pow(v1.get(i)[0] - v2.get(j)[0], 2)
                                + Math.pow(v1.get(i)[1] - v2.get(j)[1], 2));
            }
        }

        return m;

    }

    // private boolean isBorder(final Point p, final BinaryImage bi) {
    //
    // if (bi.getPixel(p.getPos(0), p.getPos(1)) == BinaryImage.OFF_VALUE) {
    // return false;
    // }
    // if (p.getPos(0) == 0 || p.getPos(0) == bi.getDimension(0) - 1 ||
    // p.getPos(1) == 0
    // || p.getPos(1) == bi.getDimension(1) - 1) {
    // return true;
    // }
    //
    // int n = 0;
    // for (int i = -1; i <= 1; i++) {
    // for (int j = -1; j <= 1; j++) {
    // if (bi.getPixel(p.getPos(0) + i, p.getPos(1) + j) ==
    // BinaryImage.OFF_VALUE) {
    // n++;
    // }
    // if (n > 1) {
    // return true;
    // }
    // }
    // }
    //
    // return false;
    //
    // }

    // public String toString() {
    // String res = "";
    // for (int i = 0; i < m_id; i++) {
    // res += m_img + " #" + i + " " + m_fn[i] + " " + m_fp[i] + " "
    // + m_total[i] + "\n";
    // }
    //
    // return res;
    //
    // }

    public Img<ByteType> getOverlayImage() {
        return m_overlayImg;
    }

}
