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
package org.knime.knip.suise.data.feat;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgView;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.iterator.IntervalIterator;
import net.imglib2.ops.operation.SubsetOperations;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

/**
 * Combines {@link PixFeatureSet}'s to one set and adds the ability to create a
 * feature image, possibly faster (if the added feature sets implement the
 * {@link FeatImgGenerator} interface.
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class PixFeatureFactory<T extends RealType<T>> {

    private final List<PixFeatureSet<T>> m_featureSetList =
            new ArrayList<PixFeatureSet<T>>();

    private final List<Integer> m_featureSetIdOffset = new ArrayList<Integer>();

    private final Collection<PixFeatureSet<T>> m_fsets;

    private int m_numFeatures = -1;

    /**
     * Creates a new feature factory
     * 
     * @param fsets
     */
    public PixFeatureFactory(final Collection<PixFeatureSet<T>> fsets) {
        m_fsets = fsets;
        updateTemporaryDataStructures();
    }

    /**
     * Updates the image on which the features are calculated.
     * 
     * @param img
     */
    void updateImg(Img<T> img) {
        for (PixFeatureSet<T> set : m_fsets) {
            set.updateImg(img);
        }
    }

    /**
     * Updates the position and orientation for which the features are
     * calculated.
     * 
     * @param dart
     */
    void updateDart(Dart dart) {
        for (PixFeatureSet<T> set : m_fsets) {
            set.updateDart(dart);
        }
    }

    /**
     * 
     * @param featID
     * @return the feature for the given feature ID, the feature id is assigned
     *         according to the order the feature sets were added
     */
    public double getFeatureValue(final int featID) {
        return m_featureSetList.get(featID).value(
                getFeatureSetFeatureID(featID));
    }

    /**
     * @return the feature values of all enabled features
     */
    public double[] getFeatureValues() {
        return getFeatureValues(new double[getNumFeatures()]);
    }

    /**
     * @return the feature values of all enabled features
     */
    public double[] getFeatureValues(final double[] vec) {
        for (int i = 0; i < vec.length; i++) {
            vec[i] =
                    m_featureSetList.get(i).value(
                            i - m_featureSetIdOffset.get(i));
        }
        return vec;
    }

    /**
     * The total number of features.
     * 
     * @return num enabled features
     */
    public int getNumFeatures() {
        return m_numFeatures;
    }

    /**
     * @param featIdx
     * @return the feature set used to calculate the feature of the given index
     */
    protected PixFeatureSet<T> getFeatureSetForFeatureIdx(final int featIdx) {
        return m_featureSetList.get(featIdx);
    }

    /**
     * Generates a feature image
     * 
     * @param srcImg the n-dimensional source image
     * @param featCalcDimionality the dimensionality of feature calculation if
     *            2, the feature calculation will be done 2D; if 3 it will be
     *            done in 3D and the used feature sets musst be able to deal
     *            with 3D images), the dimensions beyond the specified number of
     *            dimension are iterated
     * @param numOrientations the number of quantised orientations (the number
     *            of angles in 2D), if 1, no additional orientation dimension
     *            will be added to the feature image and the feature calucation
     *            is assumed to be orientation-independent
     * @return the feature image with the same dimensions as the source image
     *         plus a dimension containing the features and a dimension containg
     *         the orientations (if numOrientations > 1)
     */
    public Img<FloatType> makeFeatureImage(Img<T> srcImg,
            int featCalcDimionality, int numOrientations) {

        Img<FloatType> featImg = null;

        // 2D pixel features
        if (featCalcDimionality == 2) {

            // image dimensions + features, angles
            long[] dims =
                    new long[srcImg.numDimensions()
                            + (numOrientations == 1 ? 1 : 2)];
            for (int i = 0; i < srcImg.numDimensions(); i++) {
                dims[i] = srcImg.dimension(i);
            }
            dims[srcImg.numDimensions()] = getNumFeatures();
            if (numOrientations > 1) {
                dims[srcImg.numDimensions() + 1] = numOrientations;
            }

            featImg =
                    new ArrayImgFactory<FloatType>().create(dims,
                            new FloatType());

            // create darts
            Dart[] darts = new Dart[numOrientations];
            for (int i = 0; i < numOrientations; i++) {
                darts[i] =
                        new Dart2D(new long[srcImg.numDimensions()], i,
                                numOrientations);
            }

            if (srcImg.numDimensions() == 2) {
                generateFeatImg(srcImg, featImg, darts);
            } else {
                // create the according sub-images and iterate through the
                // dimension beyond x,y in the source and feature image

                long[] srcImgMin = new long[srcImg.numDimensions()];
                long[] srcImgMax = new long[srcImg.numDimensions()];
                srcImg.min(srcImgMin);
                srcImg.max(srcImgMax);

                long[] featImgMin = new long[featImg.numDimensions()];
                long[] featImgMax = new long[featImg.numDimensions()];

                featImg.min(featImgMin);
                featImg.max(featImgMax);

                long[] iterDims = new long[srcImg.numDimensions() - 2];
                for (int i = 2; i < srcImg.numDimensions(); i++) {
                    iterDims[i - 2] = srcImg.dimension(i);
                }

                IntervalIterator ii = new IntervalIterator(iterDims);

                while (ii.hasNext()) {
                    ii.fwd();
                    for (int i = 0; i < iterDims.length; i++) {
                        srcImgMin[i + 2] = ii.getLongPosition(i);
                        srcImgMax[i + 2] = ii.getLongPosition(i);
                        featImgMin[i + 2] = ii.getLongPosition(i);
                        featImgMax[i + 2] = ii.getLongPosition(i);
                    }

                    generateFeatImg(
                            new ImgView<T>(SubsetOperations.subsetview(srcImg,
                                    new FinalInterval(srcImgMin, srcImgMax)),
                                    srcImg.factory()),
                            SubsetOperations.subsetview(featImg,
                                    new FinalInterval(featImgMin, featImgMax)),
                            darts);

                }

            }

        } else {
            throw new UnsupportedOperationException(
                    "3D pixel feature calculation not supported, yet.");
        }
        return featImg;

    }

    private void generateFeatImg(Img<T> srcImg,
            RandomAccessibleInterval<FloatType> featImg, Dart[] darts) {

        updateImg(srcImg);

        RandomAccess<FloatType> featRA = featImg.randomAccess();

        long[] featImgMin = new long[featImg.numDimensions()];
        long[] featImgMax = new long[featImg.numDimensions()];

        int featDimIdx = featImg.numDimensions() - (darts.length > 1 ? 2 : 1);
        featImg.min(featImgMin);
        featImg.max(featImgMax);
        for (int feat = 0; feat < getNumFeatures(); feat++) {

            if (getFeatureSetForFeatureIdx(feat) instanceof FeatImgGenerator) {
                // directly let the feature set write into the feature image
                PixFeatureSet<T> fSet =
                        (PixFeatureSet<T>)getFeatureSetForFeatureIdx(feat);
                featImgMin[featDimIdx] = feat;
                featImgMax[featDimIdx] = feat;
                for (int d = 0; d < darts.length; d++) {
                    if (darts.length > 1) {
                        featImgMin[featImg.numDimensions() - 1] = d;
                        featImgMax[featImg.numDimensions() - 1] = d;
                    }
                    RandomAccessibleInterval<FloatType> subFeatImg =
                            SubsetOperations.subsetview(featImg,
                                    new FinalInterval(featImgMin, featImgMax));
                    ((FeatImgGenerator<FloatType>)fSet).generateFeatImg(
                            subFeatImg, getFeatureSetFeatureID(feat), darts[d]);
                }
            } else if (getFeatureSetForFeatureIdx(feat) instanceof WholeFeatImgGenerator) {
                // directly let the feature set write into the whole feature
                // image
                PixFeatureSet<T> fSet =
                        (PixFeatureSet<T>)getFeatureSetForFeatureIdx(feat);
                featImgMin[featDimIdx] = getFeatureSetFeatureID(feat);
                featImgMax[featDimIdx] =
                        getFeatureSetFeatureID(feat + fSet.numFeatures() - 1);
                if (darts.length > 1) {
                    featImgMin[featImg.numDimensions() - 1] = 0;
                    featImgMax[featImg.numDimensions() - 1] =
                            featImg.max(featImg.numDimensions() - 1);
                }
                RandomAccessibleInterval<FloatType> subFeatImg =
                        SubsetOperations.subsetview(featImg, new FinalInterval(
                                featImgMin, featImgMax));
                ((WholeFeatImgGenerator<FloatType>)fSet).generateFeatImg(
                        subFeatImg,
                        darts.length > 1 ? featImg.numDimensions() - 1 : -1,
                        featDimIdx);
                feat += fSet.numFeatures();
            } else {

                featRA.setPosition(feat, featDimIdx);

                for (int d = 0; d < darts.length; d++) {

                    if (darts.length > 1) {
                        featRA.setPosition(d, featImg.numDimensions() - 1);
                    }

                    Cursor<T> c = srcImg.localizingCursor();
                    while (c.hasNext()) {
                        c.fwd();
                        for (int i = 0; i < srcImg.numDimensions(); i++) {
                            featRA.setPosition(c.getIntPosition(i), i);
                        }
                        c.localize(darts[d].offset());
                        updateDart(darts[d]);
                        featRA.get().setReal(getFeatureValue(feat));
                    }
                }

            }

        }
    }

    /*
     * @return the feature id in the feature set, where featID points to
     */
    private int getFeatureSetFeatureID(final int featID) {
        return featID - m_featureSetIdOffset.get(featID);
    }

    private void updateTemporaryDataStructures() {

        int currentOffset = 0;
        for (final PixFeatureSet<T> fset : m_fsets) {
            for (int i = 0; i < fset.numFeatures(); i++) {
                m_featureSetIdOffset.add(currentOffset);
                m_featureSetList.add(fset);
            }
            currentOffset += fset.numFeatures();

        }

        m_numFeatures = 0;
        for (PixFeatureSet<T> fset : m_fsets) {
            m_numFeatures += fset.numFeatures();
        }

    }

}
