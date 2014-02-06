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

import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.real.FloatType;
import weka.clusterers.SimpleKMeans;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class ContourDataFromClusterSelection extends ContourDataExtractor {

    private final int m_numClusters;

    private final double m_maxCoverage;

    private final double m_bias;

    public ContourDataFromClusterSelection(int numClusters, double maxCoverage,
            double bias) {
        m_numClusters = numClusters;
        m_maxCoverage = maxCoverage;
        m_bias = bias;

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void extractContourData(int[] translations, int[] permutation) {
        SimpleKMeans clusterer = new SimpleKMeans();
        try {

            clusterer.setNumClusters(m_numClusters);

            // cluster the data
            ArrayList<Attribute> attInfo = new ArrayList<Attribute>();
            for (int a = 0; a < contourDataGrid().numFeatures(); a++) {
                attInfo.add(new Attribute("att" + a));
            }
            Instances data =
                    new Instances("dataset", attInfo, contourDataGrid()
                            .numVectors());
            for (double[] vec : contourDataGrid()) {
                data.add(new DenseInstance(1.0, vec));
            }
            clusterer.buildClusterer(data);

            // create clustered images p(C|x)
            Img[] imgs = new Img[m_numClusters];
            int[] dims =
                    new int[]{contourDataGrid().width(),
                            contourDataGrid().totalLength()};
            Cursor<FloatType>[] cursors = new Cursor[m_numClusters];
            for (int i = 0; i < imgs.length; i++) {
                imgs[i] =
                        new ArrayImgFactory<FloatType>().create(dims,
                                new FloatType());
                cursors[i] = imgs[i].localizingCursor();
            }

            int cluster;
            for (Instance instance : data) {
                for (int i = 0; i < cursors.length; i++) {
                    cursors[i].fwd();
                }
                cluster = clusterer.clusterInstance(instance);
                cursors[cluster].get().set(1.0f);
            }

            // greedily select the best cluster combination starting with all
            // clusters together and then removing the one whose removal
            // maximises the score of the remaining clusters
            Img<FloatType> res =
                    imgs[0].factory().create(imgs[0], new FloatType());
            Cursor<FloatType> resC = res.cursor();
            while (resC.hasNext()) {
                resC.fwd();
                resC.get().set(1.0f);
            }
            Img<FloatType> tmp = res.factory().create(res, new FloatType());

            // TODO: normalize img
            // NormalizeIterableInterval<FloatType, Img<FloatType>> imgNorm =
            // new NormalizeIterableInterval<FloatType, Img<FloatType>>();
            double score = 0;
            double bestScore = -Double.MAX_VALUE;
            double globalBestScore = -Double.MAX_VALUE;
            int bestCluster = 0;

            // ShowInSameFrame showInFrame = new ShowInSameFrame();

            for (int i = 0; i < m_numClusters; i++) {
                for (int j = 0; j < m_numClusters; j++) {
                    if (imgs[j] != null) {
                        substract(res, imgs[j], tmp);
                        score = calcScore(tmp, m_bias);
                        if (score > bestScore) {
                            bestScore = score;
                            bestCluster = j;
                        }
                    }
                }
                substract(res, imgs[bestCluster], res);
                imgs[bestCluster] = null;

                // Pair<FloatType, FloatType> minmax =
                // Operations.compute(new MinMax<FloatType>(), tmp);
                // Operations.<FloatType, FloatType> map(
                // new Normalize<FloatType>(minmax.getA().getRealDouble(),
                // minmax.getB().getRealDouble(),
                // -Float.MAX_VALUE, Float.MAX_VALUE)).compute(
                // tmp, tmp);

                // showInFrame.show(tmp, 2.0);

                if (bestScore < globalBestScore) {
                    break;
                }

                globalBestScore = bestScore;
                bestScore = -Double.MAX_VALUE;

            }

            // calculate the translations (mean positions)
            resC = res.localizingCursor();
            double meanPos = 0;
            double num = 0;
            int index = 0;
            while (resC.hasNext()) {
                resC.fwd();

                meanPos += resC.get().get() * resC.getDoublePosition(0);
                num += resC.get().get();
                index++;
                if ((index % res.dimension(0)) == 0) {
                    if (num > 0) {
                        translations[(int)((index - 1) / res.dimension(0))] =
                                (int)Math.round(meanPos / num) - CENTER_COL;
                    } else {
                        // setWeight((int)((index - 1) / res.dimension(0)), 0);
                        translations[(int)((index - 1) / res.dimension(0))] = 0;
                    }
                    meanPos = 0;
                    num = 0;
                }

            }

        } catch (Exception e) {
            // TODO Auto-generated catch block
        }

    }

    private double calcScore(Img<FloatType> img, double bias) {
        Cursor<FloatType> c = img.cursor();
        double total = 0;
        double totalWeighted = 0;
        double totalRowMax = 0;

        int index = 0;
        double max = 0;
        while (c.hasNext()) {
            c.fwd();
            total += c.get().get();
            totalWeighted += contourDataGrid().weight(index++) * c.get().get();
            max = Math.max(max, c.get().get());
            if ((index % img.dimension(0)) == 0) {
                totalRowMax += max;
                max = 0;
            }
        }

        // return (totalWeighted / total)
        // * Math.pow(totalRowMax / img.dimension(1), continuityWeight);
        if (totalRowMax / img.dimension(1) < m_maxCoverage) {
            return 0;
        } else {
            return (totalWeighted / (total + bias));
        }

    }

    private void substract(Img<FloatType> srcImg, Img<FloatType> subImg,
            Img<FloatType> resImg) {
        Cursor<FloatType> srcC = srcImg.cursor();
        Cursor<FloatType> subC = subImg.cursor();
        Cursor<FloatType> resC = resImg.cursor();
        while (srcC.hasNext()) {
            srcC.fwd();
            subC.fwd();
            resC.fwd();
            resC.get().set(srcC.get().get() - subC.get().get());
        }
    }

}
