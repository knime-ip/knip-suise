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
package org.knime.knip.suise.node.boundarymodel;

import java.awt.Polygon;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.labeling.Labeling;
import net.imglib2.ops.img.UnaryConstantRightAssignment;
import net.imglib2.ops.operation.real.binary.RealCopyLeft;
import net.imglib2.outofbounds.OutOfBoundsConstantValueFactory;
import net.imglib2.roi.IterableRegionOfInterest;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.util.ConstantUtils;
import net.imglib2.view.Views;

import org.knime.knip.core.awt.ColorLabelingRenderer;
import org.knime.knip.core.awt.Real2GreyRenderer;
import org.knime.knip.core.awt.labelingcolortable.DefaultLabelingColorTable;
import org.knime.knip.core.awt.labelingcolortable.ExtendedLabelingColorTable;
import org.knime.knip.core.awt.labelingcolortable.RandomMissingColorHandler;
import org.knime.knip.core.data.algebra.ExtendedPolygon;
import org.knime.knip.core.util.PolygonTools;
import org.knime.knip.core.util.ShowInSameFrame;
import org.knime.knip.core.util.ShowInSameFrame.ImagePlaneProducer;
import org.knime.knip.suise.data.feat.Dart;
import org.knime.knip.suise.data.feat.Dart2D;
import org.knime.knip.suise.node.boundarymodel.contourdata.ContourDataClassifier;
import org.knime.knip.suise.node.boundarymodel.contourdata.ContourDataExtractor;
import org.knime.knip.suise.node.boundarymodel.contourdata.ContourDataFromClusterSelection;
import org.knime.knip.suise.node.boundarymodel.contourdata.ContourDataFromIRI;
import org.knime.knip.suise.node.boundarymodel.contourdata.ContourDataGrid;
import org.knime.knip.suise.node.boundarymodel.contourdata.ContourDataMisc;
import org.knime.knip.suise.node.boundarymodel.contourdata.IRI;
import org.knime.knip.suise.node.boundarymodel.contourdata.VectorDataListImpl;
import org.knime.knip.suise.node.boundarymodel.contourdata.WekaContourDataClassifier;
import org.knime.knip.suise.node.boundarymodel.contourdata.WekaMIContourDataClassifier;
import org.knime.knip.suise.node.pixfeat2d.angledep.BufferedPixelFeatureSet;


/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.Classifier;
import weka.classifiers.trees.RandomForest;
import weka.clusterers.SimpleKMeans;
import weka.core.Utils;
import weka.core.matrix.Matrix;

/**
 * 
 * The data representing possible boundaries.
 * 
 * @author hornm, University of Konstanz
 */
public class BoundaryModel<T extends RealType<T>, F extends RealType<F>>
        implements Externalizable {

    /**
     * Standard deviation for the permutohedral lattice of the iterative
     * inference selection
     */
    public static final String OPTIONAL_PARAMETER_STDEV = "stdev";

    /**
     * Bias parameter needed in the cluster- and interval rule- selection
     * approaches
     */
    public static final String OPTIONAL_PARAMETER_BIAS = "bias";

    /*
     * how many times more random samples should be added to the random data
     */
    private static final int FACTOR_RANDOM_SAMPLES = 3;

    private static final int MAX_COUNTOUR_SAMPLE_POINTS = 1000;

    /*
     * the maximal allowed number of contour models for the classifcation (the
     * are chosen according to the highest score). 0 - if unlimited
     */
    // private static final int MAX_NUM_CONTOUR_MODELS = 5;

    /*
     * the number of pixels along the normalvectors taken into account
     */
    private static final int RADIUS = 10;

    private BufferedPixelFeatureSet<F> m_featSet;

    /*
     * the contour data grid the for different universes (feature sets) to
     * evaluate the model
     */
    private ContourDataGrid m_contourDataGrid = new ContourDataGrid(); // raw
                                                                       // contour
                                                                       // data

    private VectorDataListImpl m_bgData = new VectorDataListImpl(); // for each
                                                                    // universe

    /*
     * the learned classifiers for each contour cluster/model
     */
    private ContourDataClassifier m_classifier = null;

    /*
     * the resulting contour data
     */
    private ContourDataExtractor m_contourData;

    /*
     * the selected clusters for each contour and the score
     */
    private double[][] m_contourModels; // sample x cluster

    /*
     * the used weka classifier, if needed (depends on the boundary model
     * creation strategy)
     */
    private Classifier m_wekaClassifier;;

    // private Classifier m_wekaClassifier = new NaiveBayes();

    // for debugging purposes
    private ImagePlaneProducer<ByteType> m_debugImageProducer;

    private Img<T> m_contourImg;

    private static final int BG_CLUSTER_IDX_FOR_DEBUG = -1;

    private SelectionStrategy m_selectionStrategy =
            SelectionStrategy.ITERATIVE_INFERENCE_SELECTION;

    private Map<String, String> m_parameters;

    public BoundaryModel() {
        RandomForest r = new RandomForest();
        // r.setNumExecutionSlots(Runtime.getRuntime().availableProcessors());
        m_wekaClassifier = r;
        // m_wekaClassifier = new J48();
    }

    /**
     * Sets the feature factory on which basis the boundary model is build. Need
     * to be set before the samples can be added and the model build.
     * 
     * @param featFac
     */
    public void setFeatureSet(BufferedPixelFeatureSet<F> featSet) {
        m_featSet = featSet;

    }

    /**
     * Replaces the weka classifier to be used to classify images pixel-wise.
     * Depending on the boundary model, replacing it may not have no effect,
     * because it is not used
     * 
     * @param classifier
     */
    public void setWekaClassifier(Classifier classifier) {
        m_wekaClassifier = classifier;
    }

    /**
     * Sets the instance selection strategy to be used. By default
     * {@link SelectionStrategy#ITERATIVE_INFERENCE_SELECTION} will be used.
     */
    public void setSelectionStrategy(SelectionStrategy selectionStrategy) {
        m_selectionStrategy = selectionStrategy;

    }

    /**
     * Sets additional optional parameters. TODO: The available parameters and
     * their names need to be documented somewhere!!
     * 
     * @param parameters key-value-pair of parameters
     */
    public void setOptionalParameters(Map<String, String> parameters) {
        m_parameters = parameters;
    }

    /**
     * 
     * 
     * @param p the polygon representing this sample
     * @param width the source image width to randomly sample negative examples
     * @param height the source image height
     * @param srcImg the source image for debugging purposes, can be null
     * 
     */
    public void addSamples(Polygon p, int width, int height, Img<T> srcImg) {

        ExtendedPolygon extP = new ExtendedPolygon(p);

        // extP.showDebugImage((Img<T>) m_srcImg);

        extP = extP.resamplePolygon(MAX_COUNTOUR_SAMPLE_POINTS);

        m_contourData = null;
        // if (m_contourModels != null) {
        // throw new IllegalStateException(
        // "Model already build. No more samples can be added!");
        // }
        isFeatureSetSet();

        // preparing contour image for debugging if srcImg!=null
        int contourImgOffset = 0;
        if (srcImg != null && m_contourImg == null) {
            m_contourImg =
                    new ArrayImgFactory().create(
                            new int[]{RADIUS * 2 + 1, extP.length()}, srcImg
                                    .firstElement().createVariable());
        } else if (m_contourImg != null) {
            Img<T> tmp = m_contourImg.copy();
            contourImgOffset = (int)tmp.dimension(1);
            m_contourImg =
                    new ArrayImgFactory().create(new int[]{RADIUS * 2 + 1,
                            contourImgOffset + extP.length()}, srcImg
                            .firstElement().createVariable());
            Cursor<T> tmpCur = tmp.cursor();
            Cursor<T> contourImgCur = m_contourImg.cursor();
            while (tmpCur.hasNext()) {
                tmpCur.fwd();
                contourImgCur.fwd();
                contourImgCur.get().set(tmpCur.get());
            }
        }

        RandomAccess<T> contourImgRA = null;
        RandomAccess<T> srcImgRA = null;
        if (m_contourImg != null) {
            T type = m_contourImg.firstElement().createVariable();
            type.setReal(0);
            contourImgRA =
                    Views.extend(
                            m_contourImg,
                            new OutOfBoundsConstantValueFactory<T, Img<T>>(type))
                            .randomAccess();
            srcImgRA =
                    Views.extend(
                            srcImg,
                            new OutOfBoundsConstantValueFactory<T, Img<T>>(type))
                            .randomAccess();
        }

        // collect data, each point on the contour is a new instance
        int univIdx = 0;

        // retrieve contour data (data of the contour itself and in the
        // closer neighbourhood of RADIUS)
        double[][][] cData = new double[extP.length()][RADIUS * 2 + 1][];
        for (int j = 0; j < extP.length(); j++) {
            int[][] line =
                    PolygonTools.getLineAt(extP.getPointAt(j),
                            extP.getNormalVecAtPoint(j), RADIUS);
            long[] pos = new long[2];
            Dart dart =
                    new Dart2D(pos, extP.getAngleAtPoint(j),
                            m_featSet.numAngles());
            // m_featFac.updateFeatureTarget(extP.getAngleAtPoint(j));
            for (int k = 0; k < line.length; k++) {

                // if (line[k][0] < 0 || line[k][0] >= featImg.dimension(0)
                // || line[k][1] < 0
                // || line[k][1] >= featImg.dimension(1)) {
                // cData[k][j] = new double[m_featFac
                // .getNumEnabledFeatures()];
                // } else {

                // TODO: pixel feature factories must be able to deal with
                // positions beyond the image dimensions!!

                pos[0] = line[k][0];
                pos[1] = line[k][1];
                m_featSet.updateDart(dart);

                // set the feature vector
                double[] vec = new double[m_featSet.numFeatures()];
                for (int i = 0; i < vec.length; i++) {
                    vec[i] = m_featSet.value(i);
                }
                cData[j][k] = vec;

                // debug
                if (m_contourImg != null && univIdx == 0) {
                    contourImgRA.setPosition(k, 0);
                    contourImgRA.setPosition(contourImgOffset + j, 1);

                    // pixel from source image
                    srcImgRA.setPosition(line[k][0], 0);
                    srcImgRA.setPosition(line[k][1], 1);
                    contourImgRA.get().set(srcImgRA.get());

                    // pixel from feature image
                    // contourImgRA.get().setReal(m_featSet.value(22));
                }

                // }
            }

        }

        m_contourDataGrid.addContourSample(cData);

        // // negative training samples
        // extract background data
        // negative training data randomly chosen from the image
        Random r = new Random();

        long[] pos = new long[2];
        double angle;
        for (int j = 0; j < cData.length * FACTOR_RANDOM_SAMPLES; j++) {
            pos[0] = r.nextInt(width);
            pos[1] = r.nextInt(height);
            angle = r.nextDouble() * 2 * Math.PI;
            m_featSet.updateDart(new Dart2D(pos, angle, m_featSet.numAngles()));
            double[] vec = new double[m_featSet.numFeatures()];
            for (int i = 0; i < vec.length; i++) {
                vec[i] = m_featSet.value(i);
            }
            m_bgData.addVector(vec, 0);
        }
    }

    public void addSamples(IterableRegionOfInterest roi, int width, int height,
            Img<T> srcImg) {

        // roi iterable interval
        IterableInterval<BitType> ii =
                roi.getIterableIntervalOverROI(ConstantUtils.constantRandomAccessible(
                        new BitType(), roi.numDimensions()));

        // create bitmask
        Img<BitType> bitmask =
                new ArrayImgFactory<BitType>().create(ii, new BitType());

        int[] offset = new int[roi.numDimensions()];
        for (int o = 0; o < offset.length; o++) {
            offset[o] = (int)Math.round(roi.realMin(o));
        }

        Cursor<BitType> cur = ii.localizingCursor();
        RandomAccess<BitType> ra = bitmask.randomAccess();
        while (cur.hasNext()) {
            cur.fwd();
            for (int i = 0; i < bitmask.numDimensions(); i++) {
                ra.setPosition(cur.getIntPosition(i) - offset[i], i);
            }
            ra.get().set(true);
        }
        addSamples(PolygonTools.extractPolygon(bitmask, offset), width, height,
                srcImg);
    }

    /**
     * Builds the model from the afore added contour samples.
     * 
     * @throws Exception
     */
    public void buildModel() throws Exception {

        isFeatureSetSet();

        // AWTImageTools.showInFrame(m_contourImg, "contour img", 2.0);

        switch (m_selectionStrategy) {
        case ITERATIVE_INFERENCE_SELECTION:
            iterativeInferenceApproach();
            break;
        case NO_SELECTION:
            unmodifiedContourData();
            break;
        case CLUSTER_SELECTION:
            clusterSelectionApproach();
            break;
        case INTERVAL_RULE_INDUCTION:
            intervalRuleInductionApproach();
            break;
        }

        // System.out.println("objective fct: " + calcObjectiveFunction());
        // System.out.println("centrality: " + m_contourData.getCentrality());
        // System.out.println("continuity: " + m_contourData.getContinuity());

        // if (m_wekaClassifier instanceof J48) {
        // double treeSize = ((J48)m_wekaClassifier).measureTreeSize();
        // double numLeaves = ((J48)m_wekaClassifier).measureNumLeaves();
        // System.out.println("tree size: " + treeSize);
        // System.out.println("num leaves: " + numLeaves);
        // }

        // TODO: merge similar models

        // TODO: combine the informations of the different universes (e.g.
        // selecting good instances out of each universe, such that the global
        // contour remains continue)

    }

    protected void unmodifiedContourData() throws Exception {
        m_contourData = new ContourDataExtractor() {
            protected void extractContourData(int[] translations,
                    int[] permutation) {
                Arrays.fill(translations, 0);

            }
        };
        SimpleKMeans clusterer = new SimpleKMeans();
        int clustersPerSample = 2;
        clusterer.setNumClusters(clustersPerSample
                * m_contourDataGrid.numSamples());
        m_classifier =
                new WekaContourDataClassifier(m_wekaClassifier, m_contourData,
                        clusterer);
        m_classifier.buildClassifier(m_contourDataGrid, m_bgData);
        m_contourModels = new double[1][m_contourDataGrid.numClusters() + 1];
        m_contourModels =
                new double[m_contourData.numSamples()][m_contourData
                        .numClusters()];
        for (int i = 0; i < m_contourData.numVectors(); i++) {
            if (m_contourData.weight(i) > 0)
                m_contourModels[m_contourData.getSampleIndex(i)][m_contourData
                        .getClusterIdx(i)] = 1.0;
        }
        removeRedundantContourModels();
    }

    protected void iterativeInferenceApproach() throws Exception {
        double stdev = 200;
        /* Conditional random field approach */
        String val;
        if ((val = m_parameters.get(OPTIONAL_PARAMETER_STDEV)) != null) {
            stdev = Double.valueOf(val);
        }
        m_contourData = new ContourDataMisc(stdev);
        // m_contourData = new ContourDataFromCRF();
        // m_contourData = new ContourDataFromCRFNaive();
        SimpleKMeans clusterer = new SimpleKMeans();
        double clustersPerSample = 1;
        clusterer.setNumClusters(Math.max(
                1,
                (int)Math.round(clustersPerSample
                        * m_contourDataGrid.numSamples())));

        m_classifier =
                new WekaContourDataClassifier(m_wekaClassifier, m_contourData,
                        clusterer);
        m_classifier.buildClassifier(m_contourDataGrid, m_bgData);
        m_contourModels =
                new double[m_contourData.numSamples()][m_contourData
                        .numClusters()];
        for (int i = 0; i < m_contourData.numVectors(); i++) {
            if (m_contourData.weight(i) > 0)
                m_contourModels[m_contourData.getSampleIndex(i)][m_contourData
                        .getClusterIdx(i)] = 1.0;
        }
        removeRedundantContourModels();

    }

    protected void intervalRuleInductionApproach() throws Exception {
        /* Interval rule induction */
        // train the iri (interval rule induction) classifier
        double bias = 100;
        String val;
        if ((val = m_parameters.get(OPTIONAL_PARAMETER_BIAS)) != null) {
            bias = Double.valueOf(val);
        }
        final IRI miClass = new IRI();
        miClass.setBias(bias);

        // set the bias according to the mean sample length
        double meanSampleLength = 0;
        for (int i = 0; i < m_contourDataGrid.numSamples(); i++) {
            meanSampleLength += m_contourDataGrid.getSampleLength(i);
        }
        meanSampleLength /= m_contourDataGrid.numSamples();
        miClass.setBias((int)Math.round(meanSampleLength / 2));

        // extract the actual contour data to create the contour models
        m_contourData = new ContourDataFromIRI(miClass);
        // m_contourData.extractContourData(m_contourDataGrid);
        //
        new WekaMIContourDataClassifier(miClass).buildClassifier(
                m_contourDataGrid, m_bgData);

        /*
         * use the extracted contour data to feed a weka classifier
         */
        SimpleKMeans clusterer = new SimpleKMeans();
        double clustersPerSample = 1;
        clusterer.setNumClusters(Math.max(
                1,
                (int)Math.round(clustersPerSample
                        * m_contourDataGrid.numSamples())));

        m_classifier =
                new WekaContourDataClassifier(m_wekaClassifier, m_contourData,
                        clusterer);
        m_classifier.buildClassifier(m_contourDataGrid, m_bgData);
        m_contourModels =
                new double[m_contourData.numSamples()][m_contourData
                        .numClusters()];
        for (int i = 0; i < m_contourData.numVectors(); i++) {
            if (m_contourData.weight(i) > 0)
                m_contourModels[m_contourData.getSampleIndex(i)][m_contourData
                        .getClusterIdx(i)] = 1.0;
        }
        removeRedundantContourModels();

        /*
         * use this, if the retrieved interval rules should be used directly for
         * classification
         */
        // // retrieve a rule distribution for each sample, summarize
        // // distributions to cell models
        // ContourDataExtractor cd = m_contourData;
        // int numSamples = cd.numSamples();
        // final int numClusters = cd.numClusters();
        // int numVectors = cd.numVectors();
        // m_contourModels = new double[numSamples][numClusters];
        // for (int i = 0; i < numVectors; i++) {
        // m_contourModels[cd.getSampleIndex(i)][cd.getClusterIdx(i)]++;
        // }
        //
        // for (int i = 0; i < numSamples; i++) {
        // System.out.println(Arrays.toString(m_contourModels[i]));
        // Utils.normalize(m_contourModels[i],
        // m_contourModels[i][Utils.maxIndex(m_contourModels[i])]);
        // }
        //
        // // create a new classifier for each contour model
        // m_classifier = new ContourDataClassifier() {
        // public double contourProbability(double[] inst) throws Exception {
        // return 0;
        // }
        //
        // @Override
        // public void buildClassifier(ContourDataGrid cData,
        // VectorDataList bgData) throws Exception {
        // //
        // }
        //
        // public double[] contourProbDistribution(double[] inst)
        // throws Exception {
        // Instance i = new DenseInstance(1.0, inst);
        // double[] distr = new double[numClusters - 1];
        // for (int j = 0; j < distr.length; j++) {
        // distr[j] = miClass.getRule(j).distributionForInstance(i)[1];
        // }
        // return distr;
        // }
        // };

    }

    protected void clusterSelectionApproach() throws Exception {
        double bias = 100;
        String val;
        if ((val = m_parameters.get(OPTIONAL_PARAMETER_BIAS)) != null) {
            bias = Double.valueOf(val);
        }
        m_contourData = new ContourDataFromClusterSelection(10, .9, bias);
        m_classifier =
                new WekaContourDataClassifier(m_wekaClassifier, m_contourData,
                        null);
        m_classifier.buildClassifier(m_contourDataGrid, m_bgData);
        m_contourModels = new double[1][m_contourDataGrid.numClusters() + 1];
        Arrays.fill(m_contourModels[0], 1);

    }

    // protected void nnApproach() {
    // /* Neuronal network */
    // MultilayerPerceptron mp = new MultilayerPerceptron();
    // mp.setNominalToBinaryFilter(false);
    // mp.setNormalizeAttributes(false);
    // m_classifiers = new ContourDataClassifier[1];
    // m_classifiers[0] = new WekaContourDataClassifier(mp);
    // m_classifiers[0].buildClassifier(m_contourDataGrid, m_bgData);
    // m_contourModels = new double[1][m_contourDataGrid.numClusters() + 1];
    // Arrays.fill(m_contourModels[0], 1);
    // }

    /**
     * Creates an array of images with the given rectangle classified for each
     * cluster.
     * 
     * @param offset
     * 
     * @param size
     * 
     * @return
     * 
     * @throws Exception
     */
    public Img<ByteType> classifyImage(int[] offset, int[] size)
            throws Exception {

        isFeatureSetSet();

        if (m_classifier == null) {
            throw new IllegalStateException(
                    "The boundary model wasn't build, yet!");
        }

        ImgFactory<ByteType> imgFac = new ArrayImgFactory<ByteType>();
        int[] resImgSize = new int[3]; // x times y times angles
        resImgSize[0] = size[0];
        resImgSize[1] = size[1];
        resImgSize[2] = m_featSet.numAngles();
        Img<ByteType> res = imgFac.create(resImgSize, new ByteType());

        new UnaryConstantRightAssignment<ByteType, ByteType, ByteType>(
                new RealCopyLeft<ByteType, ByteType, ByteType>()).compute(res,
                new ByteType(Byte.MIN_VALUE), res);

        Cursor<ByteType> resCur = res.localizingCursor();

        long[] pos = new long[2];
        double ang;
        double fraction = 2 * Math.PI / m_featSet.numAngles();

        double[] featVec = new double[m_featSet.numFeatures()];
        while (resCur.hasNext()) {
            resCur.fwd();
            pos[0] = offset[0] + resCur.getLongPosition(0);
            pos[1] = offset[1] + resCur.getLongPosition(1);
            ang = resCur.getLongPosition(2) * fraction;
            Dart d = new Dart2D(pos, ang, m_featSet.numAngles());

            for (int i = 0; i < featVec.length; i++) {
                featVec[i] = m_featSet.value(i);
            }
            double[] distr = m_classifier.contourProbDistribution(featVec);
            double prob = 0;
            for (int i = 0; i < distr.length; i++) {
                prob = Math.max(prob, distr[i]);
            }
            resCur.get().set(
                    (byte)Math.max((byte)Math.round(prob * 255.0 - 128.0),
                            resCur.get().get()));

        }

        return res;
    }

    /**
     * @param offset
     * @param mask the interval in which the pixel will be classified. Only
     *            those pixels which are set to true will be classified, the
     *            others will be set to zero.
     * @return an image for each contour model
     * @throws Exception
     */
    public Img<ByteType>[] classifyImageContourModelwise(final int[] offset,
            final Img<BitType> mask) throws Exception {

        isFeatureSetSet();

        if (m_classifier == null) {
            throw new IllegalStateException(
                    "The boundary model wasn't build, yet!");
        }

        ImgFactory<ByteType> imgFac = new ArrayImgFactory<ByteType>();
        int[] resImgSize = new int[3]; // x times y times angles
        resImgSize[0] = (int)mask.dimension(0);
        resImgSize[1] = (int)mask.dimension(1);
        resImgSize[2] = m_featSet.numAngles();

        final int planeSize = resImgSize[0] * resImgSize[1];

        Img<ByteType>[] res = new Img[getNumModels()];
        for (int i = 0; i < res.length; i++) {
            res[i] = imgFac.create(resImgSize, new ByteType());
            for (ByteType t : res[i]) {
                t.set(Byte.MIN_VALUE);
            }
        }

        final CursorPool<Cursor<ByteType>> resCurPool =
                new CursorPool<Cursor<ByteType>>();
        resCurPool.add(res[0].localizingCursor());
        for (int i = 1; i < res.length; i++) {
            resCurPool.add(res[i].cursor());
        }

        // ExecutorService pool =
        // Executors.newFixedThreadPool(Runtime.getRuntime()
        // .availableProcessors() - 2);
        ExecutorService pool =
                Executors.newFixedThreadPool(Runtime.getRuntime()
                        .availableProcessors() + 1);

        final ContourDataClassifier c = m_classifier;

        for (int a = 0; a < m_featSet.numAngles(); a++) {

            final CursorPool<Cursor<ByteType>> curPoolCopy = resCurPool.copy();
            curPoolCopy.jumpFwd(a * planeSize);
            final Cursor<BitType> maskCur = mask.cursor();

            final BufferedPixelFeatureSet<F> featSetCopy = m_featSet.copy();

            pool.execute(new Runnable() {

                @Override
                public void run() {
                    long[] pos = new long[2];

                    double[] featVec = new double[featSetCopy.numFeatures()];

                    double prob;
                    Dart dart = null;
                    for (int i = 0; i < planeSize; i++) {
                        curPoolCopy.fwd();
                        maskCur.fwd();
                        if (maskCur.get().get()) {
                            pos[0] =
                                    offset[0]
                                            + curPoolCopy.get(0)
                                                    .getLongPosition(0);
                            pos[1] =
                                    offset[1]
                                            + curPoolCopy.get(0)
                                                    .getLongPosition(1);
                            if (dart == null
                                    || curPoolCopy.get(0).getLongPosition(2) != dart
                                            .directionIndex()) {
                                dart =
                                        new Dart2D(pos, curPoolCopy.get(0)
                                                .getIntPosition(2), m_featSet
                                                .numAngles());

                            }
                            featSetCopy.updateDart(dart);
                            for (int f = 0; f < featVec.length; f++) {
                                featVec[f] = featSetCopy.value(f);
                            }

                            double[] distr = null;
                            try {
                                distr = c.contourProbDistribution(featVec);
                            } catch (Exception e) {
                                throw new RuntimeException(e);
                            }

                            for (int m = 0; m < getNumModels(); m++) {
                                prob = 0;
                                for (int d = 0; d < distr.length; d++) {
                                    prob =
                                            Math.max(prob, distr[d]
                                                    * m_contourModels[m][d]);
                                }
                                curPoolCopy
                                        .get(m)
                                        .get()
                                        .set((byte)Math
                                                .round(prob * 255.0 - 128.0));
                            }
                        }
                    }

                }
            });
            resCurPool.jumpFwd(planeSize);
        }
        pool.shutdown();
        pool.awaitTermination(1, TimeUnit.HOURS);
        return res;

    }

    /**
     * 
     * @return the number of build contour models
     */
    public int getNumModels() {
        if (m_contourModels != null) {
            return m_contourModels.length;
        } else {
            return 0;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {

        StringBuffer sb = new StringBuffer();
        // sb.append("Contour models cluster selection:\n");
        // for (int i = 0; i < m_contourModels.length; i++) {
        // sb.append(" " + m_contourModels[i].toString() + " ("
        // + m_contourSupport[i] + ")\n");
        // }

        sb.append("\nClassifier");

        // sb.append(m_posDistr[i].toString());
        sb.append(m_classifier.toString());

        return sb.toString();

    }

    public ContourDataGrid getContourDataGrid() {
        return m_contourDataGrid;
    }

    public java.awt.Image getDebugImage() {

        BufferedImage res = null;

        if (m_contourData == null && m_contourImg != null) {

            // // create image producer
            Img<ByteType> img =
                    new ArrayImgFactory<ByteType>().create(
                            new int[]{m_contourDataGrid.width(),
                                    m_contourDataGrid.totalLength()},
                            new ByteType());
            // return new GreyImgRenderer<ByteType>().render(img, 0, 1,
            // new long[2], 1);
            m_debugImageProducer =
                    new ShowInSameFrame.ImagePlaneProducer<ByteType>(img);
            m_contourDataGrid.setContourDebugImage(m_debugImageProducer,
                    m_contourImg);
            return Toolkit.getDefaultToolkit()
                    .createImage(m_debugImageProducer);
        } else if (m_contourImg != null) {
            res =
                    new BufferedImage(m_contourDataGrid.width() * 2,
                            m_contourDataGrid.totalLength(),
                            BufferedImage.TYPE_INT_RGB);
            java.awt.Graphics g = res.getGraphics();

            Img img = m_contourData.transformContourImage(m_contourImg);
            Labeling<Integer> lab =
                    m_contourData
                            .clusterDistrLabeling(BG_CLUSTER_IDX_FOR_DEBUG);
            g.drawImage(
                    new Real2GreyRenderer<T>().render(img, 0, 1, new long[2])
                            .image(), 0, 0, null);
            ColorLabelingRenderer<Integer> labRend =
                    new ColorLabelingRenderer<Integer>();
            labRend.setLabelingColorTable(new ExtendedLabelingColorTable(
                    new DefaultLabelingColorTable(),
                    new RandomMissingColorHandler()));
            labRend.setLabelMapping(lab.firstElement().getMapping());
            g.drawImage(labRend.render(lab, 0, 1, new long[2]).image(),
                    m_contourDataGrid.width(), 0, null);
        }
        return res;
    }

    private void isFeatureSetSet() {
        if (m_featSet == null) {
            throw new IllegalStateException("No feature set set!");
        }
    }

    private void removeRedundantContourModels() {
        int numNullModels = 0;
        int numClusters = m_contourModels[0].length;
        for (int i = 0; i < m_contourModels.length; i++) {
            boolean equals;
            for (int j = i + 1; j < m_contourModels.length; j++) {
                equals = true;
                for (int k = 0; k < m_contourModels[i].length; k++) {
                    if (m_contourModels[i][k] != m_contourModels[j][k]) {
                        equals = false;
                        break;
                    }

                }
                if (equals) {
                    m_contourModels[i] = null;
                    numNullModels++;
                    break;
                }
            }
        }
        double[][] tmp =
                new double[m_contourModels.length - numNullModels][numClusters];
        int j = 0;
        for (int i = 0; i < m_contourModels.length; i++) {
            if (m_contourModels[i] == null) {
                continue;
            } else {
                tmp[j] = m_contourModels[i];
                j++;
            }
        }
        m_contourModels = tmp;
        // System.out.println("Redundant models: " + numNullModels);

    }

    private double calcLMDL(double epsilon) {

        // create sample matrix
        Matrix V =
                new Matrix(m_contourData.numFeatures(),
                        m_contourData.numVectors());

        double[] vec;
        for (int i = 0; i < m_contourData.numVectors(); i++) {
            vec = m_contourData.getVector(i);
            for (int j = 0; j < vec.length; j++) {
                V.set(j, i, vec[j]);
            }
        }

        // estimate of the covariance matrix
        Matrix W = V.times(V.transpose());

        W.times(m_contourData.numFeatures()
                / (epsilon * epsilon * m_contourData.numVectors()));

        W =
                Matrix.identity(m_contourData.numFeatures(),
                        m_contourData.numFeatures()).plus(W);

        return Utils.log2(W.det())
                * (m_contourData.numFeatures() + m_contourData.numVectors())
                / 2;

    }

    private double calcObjectiveFunction() {
        double res = 0;
        for (int i = 0; i < m_contourData.numVectors(); i++) {
            if (m_contourData.weight(i) == 0) {
                continue;
            }
            for (int j = 0; j < m_contourData.numVectors(); j++) {
                if (m_contourData.weight(j) == 0) {
                    continue;
                }
                double[] vec1 = m_contourData.getVector(i).clone();
                double[] vec2 = m_contourData.getVector(j).clone();
                for (int k = 0; k < vec2.length; k++) {
                    vec1[k] = vec1[k] / 255.0;
                    vec2[k] = vec2[k] / 255.0;
                }
                double dist = 0;
                for (int k = 0; k < vec2.length; k++) {
                    dist += Math.pow(vec1[k] - vec2[k], 2);
                }
                res += Math.exp(-dist / (2 * 10));

            }
        }
        return res;
    }

    /*
     * Serialization
     */

    /**
     * {@inheritDoc}
     */
    @Override
    public void writeExternal(ObjectOutput out) throws IOException {
        // write classifiers and clusterer
        // out.writeInt(m_numUniverses);
        // for (int i = 0; i < m_numUniverses; i++) {
        // out.writeObject(m_classifiers[i]);
        // // out.writeObject(m_posDistr[i]);
        // }
        //
        // // write contour models
        // out.writeInt(getNumModels());
        // for (int i = 0; i < getNumModels(); i++) {
        // out.writeObject(m_contourModels[i]);
        // out.writeDouble(m_contourSupport[i]);
        // }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void readExternal(ObjectInput in) throws IOException,
            ClassNotFoundException {
        // m_numUniverses = in.readInt();
        // m_classifiers = new OneClassClassifier[m_numUniverses];
        // // m_posDistr = new EM[m_numUniverses];
        // for (int i = 0; i < m_numUniverses; i++) {
        // m_classifiers[i] = (OneClassClassifier) in.readObject();
        // // m_posDistr[i] = (EM) in.readObject();
        // }
        // int numModels = in.readInt();
        // m_contourModels = new BitSet[numModels];
        // m_contourSupport = new double[numModels];
        // for (int i = 0; i < numModels; i++) {
        // m_contourModels[i] = (BitSet) in.readObject();
        // m_contourSupport[i] = in.readDouble();
        // }
    }

    public enum SelectionStrategy {
        ITERATIVE_INFERENCE_SELECTION, NO_SELECTION, CLUSTER_SELECTION, INTERVAL_RULE_INDUCTION;
    }
}

/******************************************************
 ****************** Utility classes********************
 ******************************************************/

/**
 * Pool of cursors. The underlying images of the cursors should be of same
 * dimensions and storage strategy.
 */
class CursorPool<C extends Cursor> extends ArrayList<C> {

    public void fwd() {
        for (C c : this) {
            c.fwd();
        }
    }

    public void reset() {
        for (C c : this) {
            c.reset();
        }
    }

    public void jumpFwd(final long steps) {
        for (C c : this) {
            c.jumpFwd(steps);
        }
    }

    public CursorPool<C> copy() {
        CursorPool<C> res = new CursorPool<C>();
        for (int i = 0; i < this.size(); i++) {
            res.add((C)this.get(i).copy());
        }
        return res;
    }
}
