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

import net.imglib2.Pair;
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

import org.knime.knip.core.util.ShowInSameFrame;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class ContourDataFromCRFNaive extends ContourDataExtractor {

    double[] m_prior;

    double[] m_approx0;

    double[] m_approx1;

    /**
     * {@inheritDoc}
     */
    @Override
    protected void extractContourData(int[] translations, int[] permutation) {

        ShowInSameFrame showInFrame = new ShowInSameFrame();

        m_prior =
                new double[contourDataGrid().totalLength()
                        * contourDataGrid().width()];
        Arrays.fill(m_prior, 1);
        for (int i = contourDataGrid().width() / 2; i < m_prior.length; i +=
                contourDataGrid().width()) {
            m_prior[i] = 0;
        }

        m_approx1 = new double[m_prior.length];
        m_approx0 = new double[m_approx1.length];
        Arrays.fill(m_approx0, .5);
        Arrays.fill(m_approx1, .5);

        for (int i = 0; i < m_approx1.length; i++) {
            m_approx0[i] = 1 - m_approx1[i];
        }

        for (int i = 0; i < permutation.length; i++) {
            permutation[i] = i;
        }

        ArrayImg<DoubleType, DoubleArray> debugImg =
                new ArrayImgFactory().create(new long[]{
                        contourDataGrid().width(),
                        contourDataGrid().totalLength()}, new DoubleType());

        double[] new_approx1 =
                ((DoubleArray)debugImg.update(null)).getCurrentStorageArray();
        double[] new_approx0 = new double[new_approx1.length];

        int iteration = 0;
        while (true) {

            // mean field approximation
            for (int i = 0; i < m_approx0.length; i++) {
                new_approx0[i] = 0;
                new_approx1[i] = 0;
                for (int j = 0; j < m_approx0.length; j++) {
                    if (i == j) {
                        continue;
                    }
                    // update
                    double pp = 0;
                    if (m_approx1[j] > m_approx0[j]
                            && m_approx1[i] > m_approx0[i]) {
                        pp = pairwisePotential(i, j);
                        new_approx0[i] += pp * m_approx0[j];
                        new_approx1[i] += pp * m_approx1[j];
                    }

                }

                new_approx0[i] = Math.exp(-(1 - m_prior[i]) - new_approx0[i]);
                new_approx1[i] = Math.exp(-m_prior[i] - new_approx1[i]);

                // normalization
                if (Double.isInfinite(new_approx0[i])) {
                    new_approx0[i] = 0;
                }
                if (Double.isInfinite(new_approx1[i])) {
                    new_approx1[i] = 0;
                }
                new_approx0[i] =
                        new_approx0[i] / (new_approx0[i] + new_approx1[i]);
                if (Double.isNaN(new_approx0[i])) {
                    new_approx0[i] = 0;
                }
                new_approx1[i] = 1 - new_approx0[i];
            }

            Img<FloatType> tmp = null;
            try {
                tmp =
                        debugImg.factory().imgFactory(new FloatType())
                                .create(debugImg, new FloatType());
            } catch (IncompatibleTypeException e) {
                // TODO Auto-generated catch block
            }
            new ImgConvert<DoubleType, FloatType>(debugImg.firstElement()
                    .createVariable(), new FloatType(),
                    ImgConversionTypes.DIRECT).compute(debugImg, tmp);
            Pair<FloatType, FloatType> minmax =
                    Operations.compute(new MinMax<FloatType>(), tmp);
            Operations.<FloatType, FloatType> map(
                    new Normalize<FloatType>(minmax.getA().getRealDouble(),
                            minmax.getB().getRealDouble(), -Float.MAX_VALUE,
                            Float.MAX_VALUE)).compute(tmp, tmp);
            showInFrame.show(tmp, 2.0);

            if (iteration > 20) {
                break;
            }

            iteration++;

            System.arraycopy(new_approx0, 0, m_approx0, 0, new_approx0.length);
            System.arraycopy(new_approx1, 0, m_approx1, 0, new_approx1.length);

        }

        System.out.println("num iterations: " + iteration);

        // determine translations
        for (int i = 0; i < m_approx0.length; i++) {
            if (m_approx1[i] > m_approx0[i]) {
                translations[i / contourDataGrid().width()] =
                        (i % contourDataGrid().width()) - CENTER_COL;
                new_approx1[i] = 1;
            } else {
                new_approx1[i] = 0;
            }
        }

    }

    private double pairwisePotential(int i, int j) {
        // get labels
        int labeli = getLabel(i);
        int labelj = getLabel(j);

        if (labeli == 1 && labelj == 1 && !checkGPC(i, j)) {
            // check grid position constraints
            double infinity = 100;
            return infinity;

        } else if (labeli == 1 && labelj == 1) {
            return 1 - dist(contourDataGrid().getVector(i), contourDataGrid()
                    .getVector(j));
        } else {
            return 0;
        }

    }

    // check grid position constraints
    private boolean checkGPC(int i, int j) {
        int rowi = i / contourDataGrid().width();
        int rowj = j / contourDataGrid().width();

        int coli = i % contourDataGrid().width();
        int colj = j % contourDataGrid().width();

        if (rowi == rowj) {
            return false;
        } else if ((rowi - rowj) % contourDataGrid().totalLength() == 1
                && coli - colj > 1) {
            return false;
        } else {
            return true;
        }
    }

    private int getLabel(int pos) {
        if (m_approx0[pos] > m_approx1[pos]) {
            return 0;
        } else {
            return 1;
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
        for (int i = 0; i < f2.length; i++) {
            res += Math.abs(f1[i] - f2[i]);
        }
        // res *= res;
        return Math.exp(-res / (1000000));
    }

}
