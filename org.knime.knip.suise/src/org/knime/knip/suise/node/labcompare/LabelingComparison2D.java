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
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.util.ConstantUtils;

/**
 * Compares two different images represented by their segmented components (a
 * collection of {@link Segment}-Objects) using different measures and metrics.
 * 
 * 
 */

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class LabelingComparison2D<L extends Comparable<L>, T extends IntegerType<T>> {

    private int m_width;

    private int m_height;

    private Img<ByteType> m_diffImage;

    private Img<T> m_labImg1;

    private Img<T> m_labImg2;

    // private Segment[] m_seg1Segments;
    //
    // private Segment[] m_seg2Segments;

    private int[][] m_numCommonPixels; // #ref-components x #seg-components

    private int[] m_numLab1Pixels;

    private int[] m_numLab2Pixels;

    private int m_fp;

    private int m_fn;

    private int m_hits;

    private int m_splits;

    private int m_merges;

    // the rand index
    private double m_ri;

    // the jaccard index;
    private double m_ji;

    // normalised sum of distances
    // private double m_nsd;

    // hausdorff distance (average)
    private double m_avh;

    // normalized sum of distances (averaged)
    private double m_nsd;

    // global consistency error
    private double m_gce;

    // local consistency error
    private double m_lce;

    // object-level consistency error
    private double m_oce;

    // overlapping pixels [0,1]
    // private double m_overlap;

    private static int MIN_COMMON_REGION_SIZE = 400;

    // private static double HIT_RATIO = .7;

    private List<SegmentComparison2D> m_rstatistics;

    private final NativeImgLabeling<L, T> m_lab1;

    private final NativeImgLabeling<L, T> m_lab2;

    // public final static int FOREGROUND_BLACK_SPEC = 0;
    // public final static int FOREGROUND_WHITE_SPEC = 1;

    /**
     * Note that the reference and the segmented images have to have the same
     * size and same storage strategy!!
     * 
     * @param lab1 the segmentation info retrieved from the reference image
     * @param lab2 the segmentation info retrieved from the segmented image
     * @param width width of the original images
     * @param height height of the original images
     * @param name
     */

    public LabelingComparison2D(final NativeImgLabeling<L, T> lab1,
            final NativeImgLabeling<L, T> lab2) {

        // todo: add component filter e.g. considering the component size or the
        // position (e.g. to ignore components touching the image border)
        // todo lab1.iterationOrder equals lab2.iterationOrder??

        if (lab1.dimension(0) != lab2.dimension(0)
                || lab2.dimension(1) != lab1.dimension(1)) {
            throw new IllegalArgumentException(
                    "The segmentations to compare are not of the same size!");
        }

        m_lab1 = lab1;
        m_lab2 = lab2;
        m_labImg1 = lab1.getStorageImg();
        m_labImg2 = lab2.getStorageImg();
        m_width = (int)lab1.dimension(0);
        m_height = (int)lab1.dimension(1);

        m_fn = 0;
        m_fp = 0;
        m_hits = 0;

        m_splits = 0;
        m_merges = 0;

        m_rstatistics = new Vector<SegmentComparison2D>();

        analyseComponents();

        calcRegionStatistics();

        calcMartinErrorMeasures();

        calcOCE();

        // todo: implement an evaluation metric for image segmentation of
        // multiple objects, polak, and jaccard/rand-index

    }

    public int getFPs() {
        return m_fp;
    }

    public int getFNs() {
        return m_fn;
    }

    public int getHITs() {
        return m_hits;
    }

    public int getSplits() {
        return m_splits;
    }

    public int getMerges() {
        return m_merges;
    }

    public double getRandIndex() {
        return m_ri;
    }

    public double getJaccardIndex() {
        return m_ji;
    }

    /**
     * @return the averaged hausdorff distance
     */
    public double getAvHausdorffDistance() {
        return m_avh;
    }

    /**
     * @return averaged normalized sum of distances
     */
    public double getAvgNSD() {
        return m_nsd;
    }

    public double getGCE() {
        return m_gce;
    }

    public double getLCE() {
        return m_lce;
    }

    public double getOCE() {
        return m_oce;
    }

    public Img<ByteType> getDiffImg() {
        return m_diffImage;

    }

    /*
     * Creates two maps of the components of the two different images,
     * respectively. Each array coordinate according to the respective image
     * contains either the id (map key) of the specific component covering this
     * coordinate or null.
     * 
     * Afterwards the number of pixels of each single component and the number
     * of pixels which every two components have in common (m_numCommonPixels)
     * will be count.
     * 
     * By the way an overlay image for visualization purposes will be created.
     */

    private void analyseComponents() {

        // generating the overlay image
        m_diffImage =
                new ArrayImgFactory<ByteType>().create(m_labImg1,
                        new ByteType());

        Cursor<T> lab1Cur = m_labImg1.cursor();
        Cursor<T> lab2Cur = m_labImg2.cursor();
        Cursor<ByteType> diffCur = m_diffImage.cursor();

        // overlay image
        double bg = 0;
        while (lab1Cur.hasNext()) {
            lab1Cur.fwd();
            lab2Cur.fwd();
            diffCur.fwd();
            if (lab1Cur.get().getRealDouble() != bg
                    && lab2Cur.get().getRealDouble() != bg) {
                diffCur.get().set((byte)128);
            }
            // false negative pixel
            if (lab1Cur.get().getRealDouble() == bg
                    && lab2Cur.get().getRealDouble() != bg) {
                diffCur.get().set((byte)50);
            }
            // false positive pixel
            if (lab1Cur.get().getRealDouble() != bg
                    && lab2Cur.get().getRealDouble() == bg) {
                diffCur.get().set((byte)-50);
            }

        }

        // counting the pixels in common and of each single component

        // common pixels including the background as a own component (+1)
        m_numCommonPixels =
                new int[m_lab1.getLabels().size() + 1][m_lab2.getLabels()
                        .size() + 1];
        m_numLab1Pixels = new int[m_numCommonPixels.length];
        m_numLab2Pixels = new int[m_numCommonPixels[0].length];

        lab1Cur.reset();
        lab2Cur.reset();

        int l1, l2;

        while (lab1Cur.hasNext()) {
            lab1Cur.fwd();
            lab2Cur.fwd();

            l1 = lab1Cur.get().getInteger();
            l2 = lab2Cur.get().getInteger();

            m_numCommonPixels[l1][l2]++;
            m_numLab1Pixels[l1]++;
            m_numLab2Pixels[l2]++;

        }
    }

    /*
     * Calculates regional statistics (region = connected component), which will
     * usually be summarized afterwards (e.g. average) for the whole image
     */

    private void calcRegionStatistics() {

        // calcs properties of overlapping components
        int max;
        int maxIndexSeg1 = 0;
        int maxIndexSeg2 = 0;

        int numSeg1Assignments;

        // segments without the background segment
        Segment[] segments1 = new Segment[m_numLab1Pixels.length - 1];
        Segment[] segments2 = new Segment[m_numLab2Pixels.length - 1];

        for (int i = 1; i < m_numLab2Pixels.length; i++) {
            max = 0;
            numSeg1Assignments = 0;

            for (int j = 1; j < m_numLab1Pixels.length; j++) {

                if (m_numCommonPixels[j][i] > max) {
                    max = m_numCommonPixels[j][i];
                    maxIndexSeg1 = j;
                    maxIndexSeg2 = i;
                }

                if (m_numCommonPixels[j][i] > MIN_COMMON_REGION_SIZE) {
                    numSeg1Assignments++;
                }

            }

            if (max > 1) {
                maxIndexSeg1--;
                maxIndexSeg2--;

                if (segments1[maxIndexSeg1] == null) {
                    segments1[maxIndexSeg1] = new Segment();
                    L label =
                            m_lab1.getMapping().listAtIndex(maxIndexSeg1 + 1)
                                    .get(0);
                    segments1[maxIndexSeg1].mask =
                            binaryMask(m_lab1
                                    .getIterableRegionOfInterest(label)
                                    .getIterableIntervalOverROI(
                                           ConstantUtils.constantRandomAccessible(
                                                    new BitType(), 2)));
                    segments1[maxIndexSeg1].offset = new long[2];
                    m_lab1.getExtents(label, segments1[maxIndexSeg1].offset,
                            null);
                }

                if (segments2[maxIndexSeg2] == null) {
                    segments2[maxIndexSeg2] = new Segment();
                    L label =
                            m_lab2.getMapping().listAtIndex(maxIndexSeg2 + 1)
                                    .get(0);
                    segments2[maxIndexSeg2].mask =
                            binaryMask(m_lab2
                                    .getIterableRegionOfInterest(label)
                                    .getIterableIntervalOverROI(
                                    		ConstantUtils.constantRandomAccessible(
                                                    new BitType(), 2)));
                    segments2[maxIndexSeg2].offset = new long[2];
                    m_lab2.getExtents(label, segments2[maxIndexSeg2].offset,
                            null);
                }

                SegmentComparison2D rs =
                        new SegmentComparison2D(segments1[maxIndexSeg1],
                                segments2[maxIndexSeg2]);
                // AWTImageTools.showInFrame(rs.getDebugImage(), "debug image",
                // 1);
                m_rstatistics.add(rs);

            }

            if (numSeg1Assignments > 1) {
                m_merges++;
            } else if (numSeg1Assignments == 0) {
                m_fp++;
            } else {
                m_hits++;
            }

        }

        int numSegAssignments;

        // ignore background (i=1, j=1)
        for (int i = 1; i < m_numLab1Pixels.length; i++) {
            numSegAssignments = 0;
            for (int j = 1; j < m_numLab2Pixels.length; j++) {
                if (m_numCommonPixels[i][j] > MIN_COMMON_REGION_SIZE) {
                    numSegAssignments++;

                }
            }
            if (numSegAssignments > 1) {
                m_splits++;
            } else if (numSegAssignments == 0) {
                m_fn++;
            }
        }

        m_avh = 0;
        m_nsd = 0;

        for (SegmentComparison2D rs : m_rstatistics) {
            // double fn = rs.getFNPixels();
            // double fp = rs.getFPPixels();
            // double total = rs.getTotalPixels();
            //
            // double overlap = 1 - (fn + fp) / total;
            // System.out.println("overlap: " + overlap);
            //
            // if (total < MIN_REGION_SIZE) {
            // continue;
            // }
            //
            // if (1.0 - (fn / total + fp / total) < HIT_RATIO) {
            // if (fn > fp) {
            // m_fn++;
            // } else {
            // m_fp++;
            // }
            // } else {
            // m_hits++;
            // }
            m_avh += rs.getHausdorffDistance();
            m_nsd += rs.getNSD();
        }

        m_avh /= m_rstatistics.size();
        m_nsd /= m_rstatistics.size();

    }

    /*
     * global and local consistency error (gce and lce)
     * 
     * SEE Martin, D. R.; Fowlkes, C.; Tal, D. & Malik, J. A Database of Human
     * Segmented Natural Images and its Application to Evaluating Segmentation
     * Algorithms AND Measuring Ecological Statistics EECS Department,
     * University of California, Berkeley, 2001 and Polak, M.; Zhang, H. & Pi,
     * M. An evaluation metric for image segmentation of multiple objects Image
     * and Vision Computing, 2009, 27, 1223 - 1227
     * 
     * 
     * @return the gce score, values between 0 and 1 (0 - no error, 1 - worst
     * segmentation)
     */

    private void calcMartinErrorMeasures() {
        m_gce = 0;
        m_lce = 0;

        // error between component ref[i] and seg[j]
        double p;

        // error between component seg[j] and ref[i]
        double q;

        double sp = 0;
        double sq = 0;

        // total area of intersection between ref and seg
        double n = 0;

        for (int i = 0; i < m_numLab1Pixels.length; i++) {

            for (int j = 0; j < m_numLab2Pixels.length; j++) {

                if (m_numLab1Pixels[i] == 0 || m_numLab2Pixels[j] == 0) {
                    continue;
                }

                p =
                        (1.0 - (double)m_numCommonPixels[i][j]
                                / (double)m_numLab1Pixels[i])
                                * m_numCommonPixels[i][j];
                q =
                        (1.0 - (double)m_numCommonPixels[i][j]
                                / (double)m_numLab2Pixels[j])
                                * m_numCommonPixels[i][j];
                m_lce += Math.min(p, q);
                n += m_numCommonPixels[i][j];
                sp += p;
                sq += q;
            }
        }
        m_lce /= n;
        m_gce = Math.min(sp, sq) / n;

    }

    /*
     * Calculates the object-level consistency error
     * 
     * SEE Measuring Ecological Statistics EECS Department, University of
     * California, Berkeley, 2001 and Polak, M.; Zhang, H. & Pi, M. An
     * evaluation metric for image segmentation of multiple objects Image and
     * Vision Computing, 2009, 27, 1223 - 1227
     */

    private void calcOCE() {

        double union;
        double intersec;
        double total = m_width * m_height;

        m_oce = 0;
        double oce2 = 0;
        double tmp;
        double W;
        double Wsum;

        for (int i = 0; i < m_numLab1Pixels.length; i++) {
            if (m_numLab1Pixels[i] == 0) {
                continue;
            }

            tmp = 0;
            Wsum = calcSumIntersectWithRef(i);
            for (int j = 0; j < m_numLab2Pixels.length; j++) {
                if (m_numLab2Pixels[j] == 0) {
                    continue;
                }

                intersec = m_numCommonPixels[i][j];
                union =
                        m_numLab1Pixels[i] + m_numLab2Pixels[j]
                                - m_numCommonPixels[i][j];
                W = m_numLab2Pixels[j] / Wsum;
                tmp = intersec / union * W;
            }
            m_oce += (1.0 - tmp) * (m_numLab1Pixels[i] / total);

        }

        for (int i = 0; i < m_numLab2Pixels.length; i++) {
            if (m_numLab2Pixels[i] == 0) {
                continue;
            }
            tmp = 0;
            Wsum = calcSumIntersectWithSeg(i);
            for (int j = 0; j < m_numLab1Pixels.length; j++) {
                if (m_numLab1Pixels[j] == 0) {
                    continue;
                }
                intersec = m_numCommonPixels[j][i];
                union =
                        m_numLab1Pixels[j] + m_numLab2Pixels[i]
                                - m_numCommonPixels[j][i];
                W = m_numLab1Pixels[j] / Wsum;
                tmp = intersec / union * W;
            }
            oce2 += (1.0 - tmp) * (m_numLab2Pixels[i] / total);

        }

        m_oce = Math.min(m_oce, oce2);

    }

    /*
     * calcs the cumulated area of all components intersecting with component j
     * from the reference image
     */

    private double calcSumIntersectWithRef(final int j) {
        double sum = 0;
        for (int i = 0; i < m_numLab2Pixels.length; i++) {
            if (m_numCommonPixels[j][i] > 0) {
                sum += m_numLab2Pixels[i];
            }
        }

        return sum;
    }

    /*
     * calcs the cumulated area of all components intersecting with component j
     * from the reference image
     */

    private double calcSumIntersectWithSeg(final int j) {
        double sum = 0;
        for (int i = 0; i < m_numLab1Pixels.length; i++) {
            if (m_numCommonPixels[i][j] > 0) {
                sum += m_numLab1Pixels[i];
            }
        }
        return sum;
    }

    /*
     * Calcs the cumulated area of all components
     */

    /*
     * Calculates statistics over the whole image
     */

    // private void calcStatistics() {
    // // rand and jaccard indices
    //
    // // holding the results comparing the pixel-pairs within the reference
    // // and segmentation image
    // boolean seg, ref;
    //
    // double a; // ref & seg
    // double b; // !ref & seg
    // double c; // ref & !seg
    // double d; // !ref & !seg
    //
    // a = b = c = d = 0;
    //
    // // compute 2d histogram
    //
    // // ......seg
    // // ......0 1
    // // ref 0 0 1
    // // ....1 2 3
    //
    // int[] hist2d = new int[4];
    //
    // for (int i = 0; i < m_ref.getWidth(); i++) {
    // for (int j = 0; j < m_ref.getHeight(); j++) {
    // if (m_ref.getPixel(i, j) == BinaryImage.ON_VALUE) {
    // hist2d[2]++;
    // hist2d[3]++;
    // } else {
    // hist2d[0]++;
    // hist2d[1]++;
    //
    // }
    //
    // if (m_seg.getPixel(i, j) == BinaryImage.ON_VALUE) {
    // hist2d[1]++;
    // hist2d[3]++;
    // } else {
    // hist2d[0]++;
    // hist2d[2]++;
    // }
    // }
    // }
    //
    // for (int i = 0; i < 4; i++) {
    // a += bC(hist2d[i], 2);
    // }
    // b = bC(hist2d[0] + hist2d[2], 2) + bC(hist2d[1] + hist2d[3], 2) - a;
    // c = bC(hist2d[0] + hist2d[1], 2) + bC(hist2d[2] + hist2d[3], 2) - a;
    // d = bC(m_ref.getHeight() * m_ref.getWidth(), 2) - a - b - c;
    //
    // m_ri = (a + b) / (a + b + c + d);
    // m_ji = (a + b) / (b + c + d);
    // }
    // binomial coefficient
    // private double bC(final int n, int r) {
    // double t = 1;
    //
    // int m = n - r; // r = Math.max(r, n - r);
    // if (r < m) {
    // r = m;
    // }
    //
    // for (int i = n, j = 1; i > r; i--, j++) {
    // t = t * i / j;
    // }
    //
    // return t;
    // }
    // /* Breadth-first flood fill for byte images */
    // private BinaryImage[] floodFillRegion(final int x, final int y) {
    // LinkedList<Point> q = new LinkedList<Point>();
    // q.addFirst(new Point(x, y));
    // int pixelValue;
    // int width = m_overlayImage.getWidth();
    // int height = m_overlayImage.getHeight();
    // BinaryImage tmp1 = new BinaryImage(m_ref.getWidth(),
    // m_ref.getHeight(),
    // false);
    // BinaryImage tmp2 = new BinaryImage(m_seg.getWidth(),
    // m_seg.getHeight(),
    // false);
    // Rectangle rect = new Rectangle(x, y, x + 1, y + 1);
    // while (!q.isEmpty()) {
    // Point p = q.removeLast();
    // pixelValue = m_overlayImage.getPixel(p.x, p.y);
    // if ((p.x >= 0) && (p.x < width) && (p.y >= 0) && (p.y < height)
    // && pixelValue > 0
    // && m_visited.getPixel(p.x, p.y) == BinaryImage.OFF_VALUE) {
    //
    // if (pixelValue == 100) {
    // tmp1.setPixel(p.x, p.y, BinaryImage.ON_VALUE);
    // } else if (pixelValue == 200) {
    // tmp2.setPixel(p.x, p.y, BinaryImage.ON_VALUE);
    // } else {
    // tmp1.setPixel(p.x, p.y, BinaryImage.ON_VALUE);
    // tmp2.setPixel(p.x, p.y, BinaryImage.ON_VALUE);
    // }
    //
    // rect.x = Math.min(p.x, rect.x);
    // rect.y = Math.min(p.y, rect.y);
    // rect.width = Math.max(p.x, rect.width);
    // rect.height = Math.max(p.y, rect.height);
    //
    // // label position as visited
    // m_visited.setPixel(p.x, p.y, BinaryImage.ON_VALUE);
    // if (p.x + 1 < width) {
    // q.addFirst(new Point(p.x + 1, p.y));
    // }
    // if (p.y + 1 < height) {
    // q.addFirst(new Point(p.x, p.y + 1));
    // }
    // if (p.y - 1 > 0) {
    // q.addFirst(new Point(p.x, p.y - 1));
    // }
    // if (p.x - 1 > 0) {
    // q.addFirst(new Point(p.x - 1, p.y));
    // }
    // }
    // }
    // // if (region.getPixelSet().isEmpty()) {
    // // return null;
    // // }
    // // return region;
    //
    // rect.width = rect.width - rect.x;
    // rect.height = rect.height - rect.y;
    //
    // return new BinaryImage[] { tmp1.cropRectangle(rect, null),
    // tmp2.cropRectangle(rect, null) };
    //
    // }

    private Img<BitType> binaryMask(IterableInterval<BitType> ii) {
        Img<BitType> binaryMask =
                new ArrayImgFactory<BitType>().create(ii, new BitType());
        RandomAccess<BitType> maskRA = binaryMask.randomAccess();

        Cursor<BitType> cur = ii.localizingCursor();
        while (cur.hasNext()) {
            cur.fwd();
            for (int d = 0; d < cur.numDimensions(); d++) {
                maskRA.setPosition(cur.getLongPosition(d) - ii.min(d), d);
            }
            maskRA.get().set(true);

        }
        return binaryMask;

    }
}
