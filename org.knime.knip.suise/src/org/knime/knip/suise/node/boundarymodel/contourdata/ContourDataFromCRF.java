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

import java.util.ArrayList;
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
import org.knime.knip.core.awt.AWTImageTools;
import org.knime.knip.core.data.labeling.Signature;
import org.knime.knip.core.util.ShowInSameFrame;
import org.slf4j.LoggerFactory;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class ContourDataFromCRF extends ContourDataExtractor {

    private static final double STDEV = 150.0;

    // debug
    private ShowInSameFrame m_show;

    private ArrayImg<DoubleType, DoubleArray> m_debugImg;

    private double[] m_debugArray;

    /**
     * {@inheritDoc}
     */
    @Override
    protected void extractContourData(int[] translations, int[] permutation) {

        ShowInSameFrame showInFrame = new ShowInSameFrame();

        // debug
        m_show = new ShowInSameFrame();
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

        compareSelection(prior);
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
        // int[] feat = new int[1];
        // for (int f = 0; f < contourDataGrid().numFeatures(); f++) {
        // feat[0] = f;
        instanceSelection(approx, new_approx, null);
        // determine translations
        for (int i = 0; i < approx.length; i++) {
            if (approx[i] > 0) {
                translations[i / contourDataGrid().width()] =
                        (i % contourDataGrid().width()) - CENTER_COL;

            }
        }

        for (int i = 0; i < translations.length; i++) {
            if (translations[i] == 0) {
                int j = i * contourDataGrid().width() + CENTER_COL;
                if (approx[j] == 0) {
                    setWeight(i, 0);
                }
            }

        }

        // //test
        // Img<BitType> test =
        // new ArrayImgFactory<BitType>().create(new
        // int[]{contourDataGrid().width(),
        // contourDataGrid().totalLength()}, new BitType());
        // int i = 0;
        // for (BitType t : test) {
        // if (approx[i++] == 1) {
        // t.set(true);
        // }
        // }
        // AWTImageTools.showInFrame(test, "new_approx", 1.0);

        // smooth(translations);
    }

    // } //for each feature

    /* if selectedFeatures==null all features are selected */
    private void instanceSelection(double[] approx, double[] new_approx,
            int[] selectedFeatures) {

        if (selectedFeatures == null) {
            selectedFeatures = new int[contourDataGrid().numFeatures()];
            for (int i = 0; i < selectedFeatures.length; i++) {
                selectedFeatures[i] = i;
            }
        }

        PermutohedralLattice lattice =
                new PermutohedralLattice(selectedFeatures.length, 2,
                        contourDataGrid().numVectors());

        // double[] sd = standardDeviation(m_approx);
        // Arrays.fill(m_approx, 1);
        // fadeOut(m_approx);

        int iteration = 0;
        double maxScore = -Double.MAX_VALUE;
        double lastScore = -Double.MAX_VALUE;
        ArrayList<Pair<Integer, Double>> posList =
                new ArrayList<Pair<Integer, Double>>(contourDataGrid()
                        .totalLength());
        while (true) {

            Arrays.fill(new_approx, 0);

            // re-set the priors
            // for (int i = contourDataGrid().width() / 2; i < m_prior.length; i
            // += contourDataGrid().width()) {
            // m_approx[i] = 1.0;
            // }

            if (iteration > 0) {
                lattice.beginSplatValue();
            }
            double[] tmpPos = new double[selectedFeatures.length];
            // double[] tmpPos = new double[1];
            for (int i = 0; i < contourDataGrid().numVectors(); i++) {
                // homogeneous coordinate
                double[] vec = new double[]{approx[i], 1};
                if (iteration == 0) {
                    double[] pos = contourDataGrid().getVector(i);
                    for (int j = 0; j < selectedFeatures.length; j++) {
                        tmpPos[j] = pos[selectedFeatures[j]] / (STDEV);

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
            for (int i = 0; i < contourDataGrid().numVectors(); i++) {
                lattice.slice(vec);
                new_approx[i] = vec[0] / vec[1];
                // new_approx[i] = vec[0];

                // m_approx[i] = vec[0] / vec[1];
                // min = Math.min(min, m_approx[i]);
                // max = Math.max(max, m_approx[i]);

            }

            // message passing from all xi to all xj
            // for (int i = 0; i < contourDataGrid().numVectors(); i++) {
            // double[] vec = contourDataGrid().getVector(i);
            // double[] kernDens = new double[2];
            // lattice.slice(kernDens);
            // for (int j = 0; j < contourDataGrid().numVectors(); j++) {
            // if (m_approx[j] > 0) {
            // new_approx[i] += dist(vec, contourDataGrid().getVector(j)) *
            // m_approx[j];
            //
            // }
            // }
            // new_approx[i] /= kernDens[1];
            // }

            // for (int i = 0; i < m_approx.length; i++) {
            // m_approx[i] = (m_approx[i] - min) / (max - min);
            // new_approx[i] = m_approx[i];
            // }

            // for (int i = 0; i < new_approx.length; i++) {
            // new_approx[i] -= m_approx[i];
            // }

            // retrieve row-wise maxima
            // fadeOut(new_approx);

            // find the highest probable line for each sample by dynamic
            // programming; REMARK: doesn't work very well and corrupts the
            // results
            // int accumLength = 0;
            // // Arrays.fill(m_approx, 0);
            // for (int s = 0; s < contourDataGrid().numSamples(); s++) {
            // extractMaxLine(new_approx, contourDataGrid().width(), approx,
            // accumLength, contourDataGrid()
            // .getSampleLength(s), 1);
            // accumLength += contourDataGrid().getSampleLength(s);
            // }

            rowwiseMaxima(new_approx, approx);
            // showGrid(approx, true);

            // double[] param = fitFunction(approx, new FourierSeries(), 50);
            // Arrays.fill(approx, 0);
            // evalFunction(approx, new FourierSeries(), param);
            showGrid(approx, false);

            // System.arraycopy(new_approx, 0, debugArray, 0,
            // new_approx.length);
            //
            // Img<FloatType> tmp2 =
            // new ImgConvert<DoubleType,
            // FloatType>(debugImg.firstElement().createVariable(), new
            // FloatType(),
            // ImgConversionTypes.DIRECT).compute(debugImg);
            // new NormalizeIterableInterval<FloatType,
            // Img<FloatType>>().compute(tmp2, tmp2);
            // AWTImageTools.showInFrame(tmp2, "likelihoods", 1.0);

            double posAvg = 0;
            // some additional characteristics: the sum of all weights and the
            // average position of the maxima; and normalize the new
            // approximation
            double debugSum = 0;
            int count = 0;
            for (int i = 0; i < approx.length; i++) {
                if (approx[i] > 0) {
                    debugSum += new_approx[i];
                    posAvg += i % contourDataGrid().width();
                    count++;
                }
            }
            debugSum /= count;

            // set the weight of those instances which are breaking ranks (aus
            // der Reihe tanzen) to 0
            // find the highest probable line for each sample by dynamic
            // programming
            // int accumLength = 0;
            // for (int s = 0; s < contourDataGrid().numSamples(); s++) {
            // extractMaxLine(approx, contourDataGrid().width(), new_approx,
            // accumLength, contourDataGrid()
            // .getSampleLength(s), 1);
            // accumLength += contourDataGrid().getSampleLength(s);
            // }
            //
            // for (int i = 0; i < approx.length; i++) {
            // if (approx[i] == 1) {
            // if (new_approx[i] != 1) {
            // approx[i] = 0;
            // }
            // }
            //
            // }

            //
            // //disregard those sample selection whose weight is very low and
            // which prevent the positions to average to 0 (centrality)
            // posList.clear();
            // double posSum = 0;
            // for (int i = 0; i < approx.length; i++) {
            // if (approx[i] == 1) {
            // posList.add(new Pair<Integer, Double>(i, new_approx[i]));
            // posSum += (i % contourDataGrid().width());
            // }
            // }
            // Collections.sort(posList, new Comparator<Pair<Integer, Double>>()
            // {
            // public int compare(Pair<Integer, Double> o1, Pair<Integer,
            // Double> o2) {
            // return (int)Math.round(o1.b * 1000.0 - o2.b * 1000.0);
            // }
            // });
            // int numSamples = contourDataGrid().totalLength();
            // double avg = 0;
            // double new_avg = 0;
            // while ((int)Math.round(avg) != CENTER_COL) {
            // for (int i = 0; i < posList.size(); i++) {
            // Pair<Integer, Double> p = posList.get(i);
            // new_avg = ((posSum - (p.a % contourDataGrid().width())) /
            // (double)(numSamples - 1));
            // if (Math.abs(new_avg - CENTER_COL) < Math.abs(avg - CENTER_COL))
            // {
            // avg = new_avg;
            // numSamples--;
            // posSum -= (p.a % contourDataGrid().width());
            // posList.remove(p);
            // approx[p.a] = 0;
            // break;
            // }
            // }
            // }

            // System.arraycopy(approx, 0, debugArray, 0, approx.length);
            // tmp2 =
            // new ImgConvert<DoubleType,
            // FloatType>(debugImg.firstElement().createVariable(), new
            // FloatType(),
            // ImgConversionTypes.DIRECT).compute(debugImg);
            // new NormalizeIterableInterval<FloatType,
            // Img<FloatType>>().compute(tmp2, tmp2);
            // AWTImageTools.showInFrame(tmp2, "dynamic programming", 1.0);

            // compareSelection(approx);

            System.out.println("ITERATION " + iteration + ": " + debugSum);

            // prevent the selected instances of moving towards the left
            // or right, i.e. keep them in the middle
            // posAvg /= contourDataGrid().totalLength();
            // posAvg = Math.round(CENTER_COL - posAvg);
            // int mod;
            // for (int i = 0; i < approx.length; i++) {
            // if (approx[i] == 1) {
            // approx[i] = 0;
            // mod = i % contourDataGrid().width();
            // approx[i - mod + Math.max(0, Math.min(mod + (int)posAvg,
            // contourDataGrid().width() - 1))] = 1;
            // i = i + (contourDataGrid().width() - mod);
            // }
            // }

            // check terminate conditions
            maxScore = Math.max(debugSum, maxScore);
            // if (lastScore == maxScore) { //in case of oscillations
            if (iteration == 20) {
                break;
            }
            if (debugSum < lastScore) {
                maxScore = lastScore;
            } else {
                lastScore = debugSum;
            }

            iteration++;

        }

        System.out.println("num iterations: " + iteration);
    }

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
                    tmpPos[j] = vec[j] / STDEV;
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
     * to the previous maxima-positions
     */
    private int rowwiseMaxima(double[] srcGrid, double[] resGrid) {
        double max = -1;
        double globalMax = -Double.MAX_VALUE;
        int maxIndex = 0;
        int lastMaxIndex = 0;
        int changes = 0;
        int r = contourDataGrid().width() / 2;
        for (int i = 0; i < resGrid.length; i++) {
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
        }

        // TODO make that more efficient
        // normalize
        for (int i = 0; i < resGrid.length; i++) {
            if (resGrid[i] == 1) {
                resGrid[i] = srcGrid[i] / globalMax;
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

    // private void smooth(int[] signal) {
    // int[] tmp = signal.clone();
    // for (int i = 0; i < tmp.length; i++) {
    // signal[i] = (int) Math.round(.25
    // * tmp[(i - 1 + tmp.length) % tmp.length] + .5 * tmp[i]
    // + .25 * tmp[(i + 1) % tmp.length]);
    // }
    //
    // }

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
                    LoggerFactory.getLogger(Signature.class).debug(
                            "alternative backtrack: " + l);
                    // AWTImageTools.showInFrame(m_tmp,
                    // "alt. backtrack");
                    break;
                }
            }
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
