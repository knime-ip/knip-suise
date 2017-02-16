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
package org.knime.knip.suise.node.boundarymodel.contourdata;

import java.util.Arrays;

import net.imglib2.RandomAccess;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.basictypeaccess.array.DoubleArray;
import net.imglib2.ops.operation.Operations;
import net.imglib2.ops.operation.img.unary.ImgConvert;
import net.imglib2.ops.operation.img.unary.ImgConvert.ImgConversionTypes;
import net.imglib2.ops.operation.iterableinterval.unary.MinMax;
import net.imglib2.ops.operation.real.unary.Normalize;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.optimization.fitting.CurveFitter;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;
import org.knime.core.node.NodeLogger;
import org.knime.knip.core.awt.AWTImageTools;
import org.knime.knip.core.data.labeling.Signature;
import org.knime.knip.core.util.PermutationSort;
import org.knime.knip.core.util.ShowInSameFrame;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class ContourDataMisc extends ContourDataExtractor {

    // debug
    private ShowInSameFrame m_show = new ShowInSameFrame();

    private ArrayImg<DoubleType, DoubleArray> m_debugImg;

    private double[] m_debugArray;

    private final double m_stdev;

    /**
     * @param stdev standard deviation for the permutohedral lattice
     */
    public ContourDataMisc(double stdev) {
        m_stdev = stdev;

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void extractContourData(int[] translations, int[] permutation) {

        // ShowInSameFrame showInFrame = new ShowInSameFrame();

        m_debugImg =
                new ArrayImgFactory().create(new long[]{
                        contourDataGrid().width(),
                        contourDataGrid().totalLength()}, new DoubleType());
        m_debugArray =
                ((DoubleArray)m_debugImg.update(null)).getCurrentStorageArray();

        double[] prior =
                new double[contourDataGrid().totalLength()
                        * contourDataGrid().width()];
        for (int i = contourDataGrid().width() / 2; i < prior.length; i +=
                contourDataGrid().width()) {
            prior[i] = 1.0;
        }

        // compareSelection(prior);
        double[] approx = prior.clone();

        for (int i = 0; i < permutation.length; i++) {
            permutation[i] = i;
        }

        ArrayImg<DoubleType, DoubleArray> debugImg =
                new ArrayImgFactory().create(new long[]{
                        contourDataGrid().width(),
                        contourDataGrid().totalLength()}, new DoubleType());

        double[] new_approx =
                ((DoubleArray)debugImg.update(null)).getCurrentStorageArray();

        // compareSelection(approx);

        /*
         * perform instance selection process
         */
        int from = 0;
        int to = 0;
        for (int s = 0; s < contourDataGrid().numSamples(); s++) {
            to = from + contourDataGrid().getSampleLength(s);
            instanceSelectionSimple(approx, new_approx, null, from, to);
            from = to;
        }
        // instanceSelectionSimple(approx, new_approx, null, 0,
        // contourDataGrid()
        // .totalLength() - 1);

        // determine translations
        for (int i = 0; i < approx.length; i++) {
            if (approx[i] > 0) {
                translations[i / contourDataGrid().width()] =
                        (i % contourDataGrid().width()) - CENTER_COL;

            }
        }

        // set weights of un-translated rows that don't contain any 1 to 0
        for (int i = 0; i < translations.length; i++) {
            if (translations[i] == 0) {
                int j = i * contourDataGrid().width() + CENTER_COL;
                if (approx[j] == 0) {
                    setWeight(i, 0);
                }
            }

        }
    }

    private void instanceSelectionSimple(double[] approx, double[] new_approx,
            int[] selectedFeatures, int from, int to) {

        if (selectedFeatures == null) {
            selectedFeatures = new int[contourDataGrid().numFeatures()];
            for (int i = 0; i < selectedFeatures.length; i++) {
                selectedFeatures[i] = i;
            }
        }

        int gridFrom = from * contourDataGrid().width();
        int gridTo = to * contourDataGrid().width();

        int numVectors = gridTo - gridFrom;
        PermutohedralLattice lattice =
                new PermutohedralLattice(selectedFeatures.length, 2, numVectors);

        double maxScore = -Double.MAX_VALUE;
        double lastScore = -Double.MAX_VALUE;
        int iteration = 0;
        while (true) {
            for (int i = from; i < to; i++) {
                new_approx[i] = 0;
            }

            gaussTransformHomogene(lattice, selectedFeatures, approx,
                    new_approx, iteration == 0, gridFrom, gridTo);

            rowwiseMaxima(new_approx, approx, false, gridFrom, gridTo);

            /*
             * find the highest probable line for each sample by dynamic
             * programming.
             */
            if (from == 0 && to == contourDataGrid().totalLength() - 1) {
                int accumLength = 0;
                for (int s = 0; s < contourDataGrid().numSamples(); s++) {
                    extractMaxLine(approx, contourDataGrid().width(),
                            new_approx, accumLength, contourDataGrid()
                                    .getSampleLength(s), 1);
                    accumLength += contourDataGrid().getSampleLength(s);
                }
            } else {
                extractMaxLine(approx, contourDataGrid().width(), new_approx,
                        from, to - from, 1);
            }

            /*
             * Set the weight of those instances which are breaking ranks (aus
             * der Reihe tanzen) to 0
             */
            for (int i = gridFrom; i < gridTo; i++) {
                if (approx[i] == 1) {
                    if (new_approx[i] != 1) {
                        approx[i] = 0;
                    }
                }

            }

            /*
             * or alternatively fit a function ....
             */
            // FourierSeries fct = new FourierSeries();
            // Arrays.fill(new_approx, 0);
            // double[] param = fitFunction(approx, fct, 10);
            // evalFunction(new_approx, fct, param);
            //
            // System.arraycopy(new_approx, gridFrom, approx, gridFrom, gridTo
            // - gridFrom);
            //
            // showGrid(approx, true);
            // showGrid(new_approx, false);

            // double posAvg = 0;
            // some additional characteristics: the sum of all weights and the
            // average position of the maxima; and normalize the new
            // approximation
            double debugSum = 0;
            int count = 0;
            for (int i = gridFrom; i < gridTo; i++) {
                if (approx[i] > 0) {
                    debugSum += new_approx[i];
                    // posAvg += i % contourDataGrid().width();
                    count++;
                }
            }
            debugSum /= count;

            // TODO: centralize the contour

            // check terminate conditions
            maxScore = Math.max(debugSum, maxScore);
            if (lastScore == maxScore) { // in case of oscillations
                // if (iteration == 20) {
                break;
            }
            if (debugSum < lastScore) {
                maxScore = lastScore;
            } else {
                lastScore = debugSum;
            }
            iteration++;
        }

        // System.out.println("iterations: " + iteration);

    }

    // } //for each feature

    /* if selectedFeatures==null all features are selected */
    protected void instanceSelection(double[] approx, double[] new_approx,
            int[] selectedFeatures) {

        if (selectedFeatures == null) {
            selectedFeatures = new int[contourDataGrid().numFeatures()];
            for (int i = 0; i < selectedFeatures.length; i++) {
                selectedFeatures[i] = i;
            }
        }

        PermutohedralLattice lattice =
                new PermutohedralLattice(selectedFeatures.length, 1,
                        contourDataGrid().numVectors());
        double[] homogeneCoord = new double[approx.length];
        Arrays.fill(homogeneCoord, 1.0);
        gaussTransform(lattice, selectedFeatures, homogeneCoord, homogeneCoord,
                true, true);

        // double[] sd = standardDeviation(m_approx);
        // Arrays.fill(m_approx, 1);
        // fadeOut(m_approx);

        double[] scores = new double[contourDataGrid().totalLength()];
        double maxScore = -Double.MAX_VALUE;

        for (int iteration = 0; iteration < contourDataGrid().totalLength(); iteration++) {

            Arrays.fill(new_approx, 0);
            Arrays.fill(approx, 0);
            approx[iteration * contourDataGrid().width() + CENTER_COL] = 1;

            // re-set the priors
            // for (int i = contourDataGrid().width() / 2; i < m_prior.length; i
            // += contourDataGrid().width()) {
            // m_approx[i] = 1.0;
            // }

            gaussTransform(lattice, selectedFeatures, approx, new_approx,
                    false, false);
            for (int i = 0; i < new_approx.length; i++) {
                new_approx[i] /= homogeneCoord[i];
            }

            scores[iteration] = score(new_approx);
            // rowwiseMaxima(new_approx, approx, true);

            showGrid(new_approx, false);

            if (scores[iteration] > maxScore) {
                maxScore = scores[iteration];
            }

        }

        int[] perm = PermutationSort.sort(scores);
        Arrays.fill(approx, 0);
        for (int i = 0; i < 200; i++) {
            approx[perm[perm.length - 1 - i] * contourDataGrid().width()
                    + CENTER_COL] = 1;
        }

        gaussTransform(lattice, selectedFeatures, approx, new_approx, false,
                true);
        for (int i = 0; i < new_approx.length; i++) {
            new_approx[i] /= homogeneCoord[i];
        }
        rowwiseMaxima(new_approx, approx, true, 0, contourDataGrid()
                .totalLength());

        showGrid(approx, true);

        System.out.println("test");
    }

    /*
     * Gauss tranform with the help of an permutohedral lattice for speed-up
     * WITHOUT the homogene coordinate
     */
    private void gaussTransform(PermutohedralLattice lattice,
            int[] selectedFeatures, double[] srcGrid, double[] resGrid,
            boolean splatPositions, boolean blur) {
        if (!splatPositions) {
            lattice.beginSplatValue();
        }
        double[] tmpPos = new double[selectedFeatures.length];
        // double[] tmpPos = new double[1];
        for (int i = 0; i < contourDataGrid().numVectors(); i++) {
            // homogeneous coordinate
            double[] vec = new double[]{srcGrid[i]};
            if (splatPositions) {
                double[] pos = contourDataGrid().getVector(i);
                for (int j = 0; j < selectedFeatures.length; j++) {
                    tmpPos[j] = pos[selectedFeatures[j]] / (m_stdev);

                }
                lattice.splat(tmpPos, vec);
            } else {
                lattice.splatValue(vec);
            }
        }
        // fast approx. message passing from all xi to all xj
        if (blur) {
            lattice.blur();
        }
        lattice.beginSlice();
        // double min = Double.MAX_VALUE;
        // double max = -Double.MAX_VALUE;
        double[] vec = new double[1];
        for (int i = 0; i < contourDataGrid().numVectors(); i++) {
            lattice.slice(vec);
            resGrid[i] = vec[0];
            // new_approx[i] = vec[0];

            // m_approx[i] = vec[0] / vec[1];
            // min = Math.min(min, m_approx[i]);
            // max = Math.max(max, m_approx[i]);

        }
    }

    /*
     * Gauss tranform with the help of an permutohedral lattice for speed-up
     * WITH the homogene coordinate
     */
    private void gaussTransformHomogene(PermutohedralLattice lattice,
            int[] selectedFeatures, double[] srcGrid, double[] resGrid,
            boolean splatPositions, int from, int to) {
        if (!splatPositions) {
            lattice.beginSplatValue();
        }
        double[] tmpPos = new double[selectedFeatures.length];
        // double[] tmpPos = new double[1];
        for (int i = from; i < to; i++) {
            // homogeneous coordinate
            double[] vec = new double[]{srcGrid[i], 1};
            if (splatPositions) {
                double[] pos = contourDataGrid().getVector(i);
                for (int j = 0; j < selectedFeatures.length; j++) {
                    tmpPos[j] = pos[selectedFeatures[j]] / (m_stdev);

                }
                lattice.splat(tmpPos, vec);
            } else {

                lattice.splatValue(vec);
            }
        }
        // fast approx. message passing from all xi to all xj
        lattice.blur();
        lattice.beginSlice();
        // double min = Double.MAX_VALUE;
        // double max = -Double.MAX_VALUE;
        double[] vec = new double[2];
        for (int i = from; i < to; i++) {
            lattice.slice(vec);
            resGrid[i] = vec[0] / vec[1];
            // new_approx[i] = vec[0];

            // m_approx[i] = vec[0] / vec[1];
            // min = Math.min(min, m_approx[i]);
            // max = Math.max(max, m_approx[i]);

        }
    }

    private double score(double[] grid) {

        double sum = 0;
        double wSum = 0;
        double bias = 50;
        int j;
        for (int i = 0; i < grid.length; i++) {
            j = Math.abs(i % contourDataGrid().width() - CENTER_COL);
            sum += grid[i];
            wSum += grid[i] * ((double)(CENTER_COL - j) / CENTER_COL);
        }

        return wSum / (sum + bias);

        // double weightedSum = 0;
        // double sum = 0;
        // for (int i = 0; i < grid.length; i++) {
        // weightedSum += grid[i] * (i % contourDataGrid().width() -
        // CENTER_COL);
        // sum += grid[i];
        // }
        // return Math.abs(weightedSum / sum);

    }

    /*
     * Shows the grid in a window
     */
    private void showGrid(double[] grid, boolean newWindow) {
        System.arraycopy(grid, 0, m_debugArray, 0, grid.length);

        Img<FloatType> tmp = null;
        try {
            tmp =
                    m_debugImg.factory().imgFactory(new FloatType())
                            .create(m_debugImg, new FloatType());
        } catch (IncompatibleTypeException e) {
            // TODO Auto-generated catch block
        }
        new ImgConvert<DoubleType, FloatType>(m_debugImg.firstElement()
                .createVariable(), new FloatType(), ImgConversionTypes.DIRECT)
                .compute(m_debugImg, tmp);

        Pair<FloatType, FloatType> minmax =
                Operations.compute(new MinMax<FloatType>(), tmp);
        Operations.<FloatType, FloatType> map(
                new Normalize<FloatType>(minmax.getA().getRealDouble(), minmax
                        .getB().getRealDouble(), -Float.MAX_VALUE,
                        Float.MAX_VALUE)).compute(tmp, tmp);
        if (newWindow) {
            AWTImageTools.showInFrame(tmp, "debug");
        } else {
            m_show.show(tmp, 1.0);
        }

    }

    /*
     * Fits a parametric function to the grid
     */
    private double[] fitFunction(double[] approx,
            ParametricUnivariateFunction function, int numParam) {
        CurveFitter fitter = new CurveFitter(new LevenbergMarquardtOptimizer());
        // CurveFitter fitter = new CurveFitter(new GaussNewtonOptimizer());
        double x;
        double y;
        for (int i = 0; i < approx.length; i++) {
            if (approx[i] > 0) {
                x = i / contourDataGrid().width();
                x = x / contourDataGrid().totalLength() * 2 * Math.PI - Math.PI;
                y = i % contourDataGrid().width() - CENTER_COL;
                fitter.addObservedPoint(approx[i], x, y);
            }

        }

        return fitter.fit(function, new double[numParam]);

    }

    /*
     * Evaluates a parametric function and draws it on the grid
     */
    private void evalFunction(double[] approx,
            ParametricUnivariateFunction function, double[] param) {
        double x;
        double y;
        for (int i = 0; i < contourDataGrid().totalLength(); i++) {
            x =
                    (i / (double)contourDataGrid().totalLength()) * 2 * Math.PI
                            - Math.PI;
            y = function.value(x, param);
            if (Math.abs(y) >= CENTER_COL) {
                continue;
            }
            approx[i * contourDataGrid().width()
                    + ((int)Math.round(y + CENTER_COL))] = 1;
        }

    }

    private ShowInSameFrame compareSelection = new ShowInSameFrame();

    private void compareSelection(double[] grid) {
        PermutohedralLattice lat =
                new PermutohedralLattice(contourDataGrid().numFeatures(),
                        contourDataGrid().totalLength() + 1, contourDataGrid()
                                .totalLength());
        double[] value = new double[contourDataGrid().totalLength() + 1];
        double[] tmpPos = new double[contourDataGrid().numFeatures()];
        int i = 0;
        int numPos = 0;
        for (int g = 0; g < grid.length; g++) {
            if (grid[g] == 1) {
                Arrays.fill(value, 0);
                value[i] = 1;
                value[value.length - 1] = 1;
                double[] vec = contourDataGrid().getVector(g);
                for (int j = 0; j < vec.length; j++) {
                    tmpPos[j] = vec[j] / m_stdev;
                }
                lat.splat(tmpPos, value);
                i++;
                numPos++;
            }

        }
        lat.blur();
        lat.beginSlice();
        Img<FloatType> img =
                new ArrayImgFactory<FloatType>().create(new int[]{value.length,
                        value.length}, new FloatType());
        RandomAccess<FloatType> tmpRA = img.randomAccess();
        for (i = 0; i < numPos; i++) {
            lat.slice(value);
            tmpRA.setPosition(i, 1);

            for (int j = 0; j < value.length - 1; j++) {
                tmpRA.setPosition(j, 0);
                // double val = value[j] / value[value.length - 1];
                double val = value[j];
                // value[value.length - 1] += val;
                tmpRA.get().setReal(val);
            }
            value[value.length - 1] = 0;
            // tmpRA.setPosition(value.length - 1, 0);
            // tmpRA.get().setReal(value[value.length - 1]);

        }
        Pair<FloatType, FloatType> minmax =
                Operations.compute(new MinMax<FloatType>(), img);
        Operations.<FloatType, FloatType> map(
                new Normalize<FloatType>(minmax.getA().getRealDouble(), minmax
                        .getB().getRealDouble(), -Float.MAX_VALUE,
                        Float.MAX_VALUE)).compute(img, img);
        compareSelection.show(img, 1.0);
    }

    /*
     * retrieves the row-wise maxima und returns the number of changes compared
     * to the previous maxima-positions, if mask is set to true, the maxima will
     * be determined but the original value taken at the position instead of "1"
     */
    private int rowwiseMaxima(double[] srcGrid, double[] resGrid, boolean mask,
            int from, int to) {
        double max = -1;
        double globalMax = -Double.MAX_VALUE;
        double globalMin = Double.MAX_VALUE;
        int maxIndex = from;
        int lastMaxIndex = from;
        int changes = 0;
        int r = contourDataGrid().width() / 2;
        for (int i = from; i < to; i++) {
            if (i % contourDataGrid().width() == 0) {
                max = -1;
                resGrid[maxIndex] = 1;
                if (maxIndex != lastMaxIndex) {
                    changes++;
                }

            }
            if (srcGrid[i] > max) {
                max = srcGrid[i];
                maxIndex = i;
                globalMax = Math.max(max, globalMax);
            }
            // else if (srcGrid[i] == max
            // && Math.abs(i % contourDataGrid().width() - r) <
            // Math.abs(maxIndex % contourDataGrid().width() - r)) {
            // max = srcGrid[i];
            // maxIndex = i;
            // }

            if (resGrid[i] > 0) {
                resGrid[i] = 0;
                lastMaxIndex = i;
            }
            globalMin = Math.min(globalMin, srcGrid[i]);
        }

        // TODO make that more efficient
        // normalize
        if (mask) {
            for (int i = from; i < to; i++) {
                if (resGrid[i] == 1) {
                    resGrid[i] =
                            (srcGrid[i] - globalMin) / (globalMax - globalMin);
                }
            }
        }

        return changes;
    }

    private double[] standardDeviation(double[] grid) {
        double[] sum = new double[contourDataGrid().numFeatures()];
        double[] sumSquare = new double[contourDataGrid().numFeatures()];
        double count = 0;
        for (int i = 0; i < grid.length; i++) {
            if (grid[i] == 1) {
                double[] vec = contourDataGrid().getVector(i);
                for (int j = 0; j < sum.length; j++) {
                    sum[j] += vec[j];
                    sumSquare[j] += vec[j] * vec[j];
                    count++;
                }
            }
        }

        for (int i = 0; i < sum.length; i++) {
            sum[i] =
                    Math.sqrt((sumSquare[i] - ((sum[i] * sum[i]) / count))
                            / (count - 1));
        }
        return sum;
    }

    private void normalizeRow(double[] srcGrid, double[] resGrid) {
        double max = -Double.MAX_VALUE;
        double min = Double.MAX_VALUE;

        boolean firstPass = false;

        for (int i = 0; i < resGrid.length; i++) {
            if (i % contourDataGrid().width() == 0) {
                if (firstPass) {
                    i -= contourDataGrid().width();
                } else {
                    max = -Double.MAX_VALUE;
                    min = Double.MAX_VALUE;
                }
                firstPass = !firstPass;

            }

            if (firstPass) {
                max = Math.max(srcGrid[i], max);
                min = Math.min(srcGrid[i], min);
            } else {
                resGrid[i] = (srcGrid[i] - min) / (max - min);
            }

        }
    }

    private void fadeOut(double[] grid) {
        double r = contourDataGrid().width() / 2;
        for (int i = 0; i < grid.length; i++) {
            double col = i % contourDataGrid().width() - r;
            grid[i] *= (r - Math.abs(col)) / r;
        }
    }

    private int getVecIdx(int rowIdx, double[] list) {
        int offset = rowIdx * contourDataGrid().width();
        for (int i = offset; i < offset + contourDataGrid().width(); i++) {
            if (list[i] == 1) {
                return i;
            }
        }
        return -1;
    }

    private double dist(double[] f1, double[] f2) {
        double res = 0;
        // int[] sel = new int[]{4, 5, 6, 7, 13, 14, 15, 16};
        int[] sel = new int[]{f1.length / 2};
        // int[] sel = new int[f2.length];
        // for (int i = 0; i < f2.length; i++) {
        // sel[i] = i;
        // }
        for (int i = 0; i < sel.length; i++) {
            res += Math.abs(f1[sel[i]] - f2[sel[i]]);
        }
        // res *= res;
        return Math.exp(-res / 10000000);
    }

    /*
     * Helper for the constructor to extract the Signature-Line by dynamic
     * programming. The polar image is considered as a graph and we are
     * searching for the "best CLOSED path" (in this case longest path). In each
     * step the best of the available parents will be choosen, his weighted
     * added to the current node and the path direction stored. ...
     * 
     * res - the result as a grid with the selected samples set to 1; posRes -
     * the result as a list of column positions
     */

    private void extractMaxLine(double[] weights, int width, double[] res,
            int from, int length, int maxLineVariance) {

        double[][] scores = new double[length][width];
        int fromIdx = from * width;
        for (int i = fromIdx; i < fromIdx + width; i++) {
            scores[0][i - fromIdx] = weights[i];
        }
        // scores[0] = weights[0].clone();

        int[][] dirs = new int[length][width];

        // calc the shortest pathes
        double max;
        int dir = 0;

        for (int i = 1; i < scores.length; i++) {
            for (int j = 0; j < scores[0].length; j++) {

                // the maximum within the max line variance
                max = scores[i - 1][j];
                dir = 0;

                for (int s = -maxLineVariance; s <= maxLineVariance; s++) {
                    if (s == 0 || j + s < 0 || j + s > scores[0].length - 1) {
                        continue;
                    }
                    if (scores[i - 1][j + s] > max) {
                        max = scores[i - 1][j + s];
                        dir = s;
                    }
                }
                scores[i][j] = max + weights[fromIdx + i * width + j];
                dirs[i][j] = dir;
            }
        }

        // backtrack the best path two times wich in the most cases
        // results in a
        // closed contour. Else, directly find the maximal closed path.
        // int[] res = new int[weights.length];
        dir = dirs[0].length / 2;
        int firstDir = 0;
        int lastDir = 0;
        for (int i = fromIdx; i < fromIdx + length * width; i++) {
            res[i] = 0;
        }
        // for (int i = from; i < from + length; i++) {
        // posRes[i] = 0;
        // }

        for (int count = 0; count < 2; count++) {

            res[fromIdx + (length - 1) * width + dir] = 1;
            firstDir = dir;
            // posRes[posRes.length - 1] = dir;
            for (int i = dirs.length - 2; i >= 0; i--) {
                dir += dirs[i + 1][dir];
                res[i * width + fromIdx + dir] = 1;
                lastDir = dir;
                // posRes[from + i] = dir;
            }
            // if (Math.abs(res[0] - res[res.length - 1]) <= maxLineVariance) {
            if (Math.abs(lastDir - firstDir) <= maxLineVariance) {
                break;
            } else {
                for (int i = fromIdx; i < fromIdx + length * width; i++) {
                    res[i] = 0;
                }
                // for (int i = from; i < from + length; i++) {
                // posRes[i] = 0;
                // }
            }
        }
        // score = scores[scores.length - 1][res[res.length - 1]] /
        // scores.length;

        // if the found best path isn't closed, backtrack all possible
        // leaves
        // until a closed path was found

        // if (Math.abs(res[0] - res[res.length - 1]) > maxLineVariance) {
        if (Math.abs(lastDir - firstDir) > maxLineVariance) {
            // collect the leaves
            max = 0;
            IndexedDouble[] leaves = new IndexedDouble[scores[0].length];
            for (int i = 0; i < dirs[0].length; i++) {
                leaves[i] = new IndexedDouble(scores[scores.length - 1][i], i);
            }

            Arrays.sort(leaves);

            // backtrack to retrieve the path from the best leaf as
            // long as a
            // closed path was found
            for (int l = scores[0].length - 1; l >= 0; l--) {
                dir = leaves[l].getIndex();
                res[fromIdx + (length - 1) * width + dir] = 1;
                firstDir = dir;
                // posRes[dirs.length - 1] = dir;
                for (int i = dirs.length - 2; i >= 0; i--) {
                    dir += dirs[i + 1][dir];
                    // posRes[from + i] = dir;
                    res[i * width + fromIdx + dir] = 1;
                    lastDir = dir;

                }
                // score = leaves[l].getVal() / scores.length;
                // if (Math.abs(res[0] - res[res.length - 1]) <=
                // maxLineVariance) {
                if (Math.abs(lastDir - firstDir) <= maxLineVariance) {
                    NodeLogger.getLogger(Signature.class).debug(
                            "alternative backtrack: " + l);
                    // AWTImageTools.showInFrame(m_tmp,
                    // "alt. backtrack");
                    break;
                }
            }
        }
    }

    /*
     * Fourer Series function without the y-offset
     */
    private class FourierSeries implements ParametricUnivariateFunction {

        public double value(double x, double... par) {
            double y = 0;
            for (int i = 0; i < par.length; i++) {
                y += par[i] * Math.cos(i * x);
                i++;
                y += par[i] * Math.sin(i * x);
            }
            return y;
        }

        public double[] gradient(double x, double... par) {
            double[] gradient = new double[par.length];
            for (int i = 0; i < par.length; i++) {
                gradient[i] = Math.cos(i * x);
                i++;
                gradient[i] = Math.sin(i * x);
            }
            return gradient;

        }
    }

    /*
     * Hepler to associate a double value with an index. Here: to keep the
     * original index of a double list after sorting.
     */
    private static class IndexedDouble implements Comparable<IndexedDouble> {

        double m_val;

        int m_index;

        public IndexedDouble(double val, int index) {
            m_val = val;
            m_index = index;
        }

        // public double getVal() {
        // return m_val;
        // }

        public int getIndex() {
            return m_index;
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public int compareTo(IndexedDouble val2) {
            return (int)Math.round(m_val * 100 - val2.m_val * 100);
        }
    }

}
