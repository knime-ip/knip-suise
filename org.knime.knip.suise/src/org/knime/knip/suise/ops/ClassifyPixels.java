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

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.iterator.IntervalIterator;
import net.imglib2.ops.img.UnaryObjectFactory;
import net.imglib2.ops.img.UnaryOperationAssignment;
import net.imglib2.ops.operation.UnaryOutputOperation;
import net.imglib2.ops.operation.real.unary.RealConstant;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.numeric.RealType;

import org.knime.knip.core.KNIPGateway;

/**
 * TODO Auto-generated 
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
import weka.classifiers.Classifier;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;

/**
 * 
 * @author Martin Horn, University of Konstanz
 */
public class ClassifyPixels<T extends RealType<T>, RT extends RealType<RT>>
		implements UnaryOutputOperation<Img<T>, Img<RT>[]> {

	private final boolean m_createLabeling;

	private RandomAccessibleInterval<LabelingType<String>> m_labeling;

	private final Classifier m_classifier;

	private final int m_numClasses;

	private final int[] m_dimIndices;

	private final int m_featDimIdx;

	private final Instances m_dataset;

	private final RT m_resType;

	public ClassifyPixels(Classifier classifier, Instances dataset,
			int numClasses, int[] dimIndices, int featDimIdx,
			boolean createLabeling, RT resType) {
		m_classifier = classifier;
		m_dataset = dataset;
		m_numClasses = numClasses;
		m_dimIndices = dimIndices;
		m_featDimIdx = featDimIdx;
		m_createLabeling = createLabeling;
		m_resType = resType;

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Img<RT>[] compute(Img<T> op, Img<RT>[] r) {

		RandomAccess<LabelingType<String>> labelingRA = null;
		if (m_createLabeling) {
			long[] dims = new long[m_dimIndices.length];
			for (int j = 0; j < dims.length; j++) {
				dims[j] = op.dimension(j);
			}

			m_labeling = (RandomAccessibleInterval<LabelingType<String>>) KNIPGateway
					.ops().createImgLabeling(dims);
			labelingRA = m_labeling.randomAccess();
		}

		long[] max = new long[m_dimIndices.length];
		for (int d = 0; d < m_dimIndices.length; d++) {
			max[d] = op.max(m_dimIndices[d]);
		}

		IntervalIterator ii = new IntervalIterator(
				new long[m_dimIndices.length], max);

		RandomAccess<RT>[] resRAs = new RandomAccess[r.length];
		for (int i = 0; i < resRAs.length; i++) {
			resRAs[i] = r[i].randomAccess();
		}
		RandomAccess<T> imgRA = op.randomAccess();

		double[] featVec = new double[(int) op.dimension(m_featDimIdx)];
		Instance instance = new DenseInstance(1.0, featVec);
		instance.setDataset(m_dataset);
		while (ii.hasNext()) {
			ii.fwd();
			for (int d = 0; d < m_dimIndices.length; d++) {
				imgRA.setPosition(ii.getLongPosition(d), m_dimIndices[d]);
			}

			for (int f = 0; f < op.dimension(m_featDimIdx); f++) {
				imgRA.setPosition(f, m_featDimIdx);
				featVec[f] = imgRA.get().getRealDouble();
			}

			double[] probs = null;
			try {
				probs = m_classifier.distributionForInstance(instance);

				if (m_createLabeling) {
					for (int d = 0; d < m_dimIndices.length; d++) {
						labelingRA.setPosition(ii.getLongPosition(d), d);
						labelingRA.get().clear();
						labelingRA.get().add(
								m_dataset.classAttribute().value(
										Utils.maxIndex(probs)));
					}
				}

			} catch (Exception e) {
				e.printStackTrace();
			}
			for (int i = 0; i < probs.length; i++) {
				if (probs[i] > 0) {
					for (int d = 0; d < m_dimIndices.length; d++) {
						resRAs[i].setPosition(ii.getLongPosition(d), d);
						resRAs[i].get().setReal(
								probs[i]
										* (m_resType.getMaxValue() - m_resType
												.getMinValue())
										+ m_resType.getMinValue());
					}
				}
			}

		}

		return r;

	}

	public RandomAccessibleInterval<LabelingType<String>> getLabeling() {
		return m_labeling;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public UnaryOutputOperation<Img<T>, Img<RT>[]> copy() {
		return new ClassifyPixels<T, RT>(m_classifier, m_dataset, m_numClasses,
				m_dimIndices, m_featDimIdx, m_createLabeling, m_resType);
	}

	/**
	 * {@inheritDoc}
	 */
	public UnaryObjectFactory<Img<T>, Img<RT>[]> bufferFactory() {
		return new UnaryObjectFactory<Img<T>, Img<RT>[]>() {

			@Override
			public Img<RT>[] instantiate(Img<T> a) {
				Img<RT>[] res = new Img[m_numClasses];
				UnaryOperationAssignment<RT, RT> set = new UnaryOperationAssignment<RT, RT>(
						new RealConstant<RT, RT>(m_resType.getMinValue()));
				for (int i = 0; i < res.length; i++) {
					try {
						long[] dims = new long[m_dimIndices.length];
						for (int j = 0; j < dims.length; j++) {
							dims[j] = a.dimension(j);
						}
						res[i] = a.factory().imgFactory(m_resType)
								.create(dims, m_resType);

						set.compute(res[i], res[i]);
					} catch (IncompatibleTypeException e) {
						e.printStackTrace();
					}
				}
				return res;
			}
		};
	}

}
