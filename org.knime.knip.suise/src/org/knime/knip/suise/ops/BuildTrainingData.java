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
package org.knime.knip.suise.ops;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.ops.operation.BinaryObjectFactory;
import net.imglib2.ops.operation.BinaryOutputOperation;
import net.imglib2.outofbounds.OutOfBounds;
import net.imglib2.outofbounds.OutOfBoundsBorder;
import net.imglib2.roi.RectangleRegionOfInterest;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.numeric.RealType;

import org.knime.core.node.NodeLogger;

import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

/**
 * 
 * @param <L>
 * @param <T>
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class BuildTrainingData<L extends Comparable<L>, T extends RealType<T>>
		implements
		BinaryOutputOperation<RandomAccessibleInterval<LabelingType<L>>, Img<T>, Instances> {

	private final int[] m_dimIndices;

	private final int m_featDim;

	private final List<String> m_classLabels;

	private final double m_samplingRate;

	private final boolean m_balanceInstancePerClass;

	public BuildTrainingData(List<String> classLabels, int[] dimIndices,
			int featDim, double samplingRate, boolean balanceInstancePerClass) {
		m_classLabels = classLabels;
		m_dimIndices = dimIndices;
		m_featDim = featDim;
		m_samplingRate = samplingRate;
		m_balanceInstancePerClass = balanceInstancePerClass;

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Instances compute(RandomAccessibleInterval<LabelingType<L>> lab,
			Img<T> img, Instances r) {
		Random rand = new Random();

		double[] extent = new double[lab.numDimensions()];
		for (int d = 0; d < m_dimIndices.length; d++) {
			extent[m_dimIndices[d]] = lab.max(m_dimIndices[d]);
		}
		RectangleRegionOfInterest roi = new RectangleRegionOfInterest(
				new double[lab.numDimensions()], extent);

		Cursor<LabelingType<L>> labCur = roi.getIterableIntervalOverROI(lab)
				.localizingCursor();
		OutOfBounds<T> imgRA = new OutOfBoundsBorder<T>(img);

		LabelRegions<L> regions = new LabelRegions<L>(lab);
		// get the class distributions
		Map<L, Double> classDistr = null;
		if (m_balanceInstancePerClass) {
			long sum = 0;
			long area;
			Collection<L> labels = regions.getExistingLabels();
			classDistr = new HashMap<L, Double>(labels.size());
			for (L label : labels) {
				area = regions.getLabelRegion(label).size();
				sum += area;
				classDistr.put(label, new Double(area));
			}
			// determine the new sampling rate for each class individually
			double instancesPerClass = (double) sum / (double) labels.size();
			for (L label : labels) {
				Double sampleRate = instancesPerClass / classDistr.get(label)
						* m_samplingRate;
				classDistr.put(label, sampleRate);
			}
		}

		long[] tmpPos = new long[imgRA.numDimensions()];
		while (labCur.hasNext()) {
			labCur.fwd();
			for (int d = 0; d < m_dimIndices.length; d++) {
				imgRA.setPosition(labCur.getLongPosition(m_dimIndices[d]),
						m_dimIndices[d]);
				if (imgRA.isOutOfBounds()) {
					imgRA.localize(tmpPos);
					NodeLogger.getLogger(getClass()).warn(
							"Labeling reaches beyond the feature image. Position "
									+ Arrays.toString(tmpPos) + " skipped.");
					continue;
				}

			}
			if (!labCur.get().isEmpty()) {

				if (m_balanceInstancePerClass) {
					if (rand.nextDouble() >= classDistr.get(labCur.get()
							.iterator().next())) {
						continue;
					}
				} else {
					if (rand.nextDouble() >= m_samplingRate) {
						continue;
					}
				}

				double[] featVec = new double[(int) img.dimension(m_featDim)];
				for (int f = 0; f < img.dimension(m_featDim); f++) {
					imgRA.setPosition(f, m_featDim);
					featVec[f] = imgRA.get().getRealDouble();
				}
				for (L classLabel : labCur.get()) {
					Instance instance = new DenseInstance(1.0, featVec);
					instance.insertAttributeAt(instance.numAttributes());
					instance.setDataset(r);
					instance.setClassValue(classLabel.toString());

					r.add(instance);

				}
			}
		}
		return r;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public BinaryOutputOperation<RandomAccessibleInterval<LabelingType<L>>, Img<T>, Instances> copy() {
		return new BuildTrainingData<L, T>(m_classLabels, m_dimIndices,
				m_featDim, m_samplingRate, m_balanceInstancePerClass);
	}

	/**
	 * {@inheritDoc}
	 */
	public BinaryObjectFactory<RandomAccessibleInterval<LabelingType<L>>, Img<T>, Instances> bufferFactory() {
		return new BinaryObjectFactory<RandomAccessibleInterval<LabelingType<L>>, Img<T>, Instances>() {

			@Override
			public Instances instantiate(
					RandomAccessibleInterval<LabelingType<L>> inputA,
					Img<T> inputB) {
				// build training set
				ArrayList<Attribute> attr = new ArrayList<Attribute>();
				for (int a = 0; a < inputB.dimension(m_featDim); a++) {
					attr.add(new Attribute("attr" + a));
				}
				Instances instances = new Instances("data", attr,
						m_classLabels.size() * 20);

				instances.insertAttributeAt(new Attribute("class",
						m_classLabels), instances.numAttributes());
				instances.setClassIndex(instances.numAttributes() - 1);
				return instances;
			}

		};
	}
}
