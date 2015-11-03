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

import net.imglib2.img.Img;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.type.numeric.integer.ByteType;

import org.knime.knip.core.features.FeatureSet;
import org.knime.knip.core.features.FeatureTargetListener;

/**
 * 
 * @author <a href="mailto:horn_martin@gmx.de">Martin Horn</a>
 */
public class LabelingCompareFeatureSet<L extends Comparable<L>> implements
        FeatureSet {

    LabelingComparison2D m_labelingComparison;

    public static final String[] FEATURES = new String[]{"Added", "Missing",
            "HITs", "Splits", "Merges", "RandIndex", "JaccardIndex",
            "Avg Hausdorff Distance", "NSD", "GCE", "LCE", "OCE"};

    /**
     * {@inheritDoc}
     */
    @Override
    public double value(int id) {
        switch (id) {
        case 0:
            return m_labelingComparison.getFPs();
        case 1:
            return m_labelingComparison.getFNs();
        case 2:
            return m_labelingComparison.getHITs();
        case 3:
            return m_labelingComparison.getSplits();
        case 4:
            return m_labelingComparison.getMerges();
        case 5:
            return m_labelingComparison.getRandIndex();
        case 6:
            return m_labelingComparison.getJaccardIndex();
        case 7:
            return m_labelingComparison.getAvHausdorffDistance();
        case 8:
            return m_labelingComparison.getAvgNSD();
        case 9:
            return m_labelingComparison.getGCE();
        case 10:
            return m_labelingComparison.getLCE();
        case 11:
            return m_labelingComparison.getOCE();
        default:
            return 0;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void enable(int id) {

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String name(int id) {
        return FEATURES[id];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int numFeatures() {
        return FEATURES.length;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String featureSetId() {
        return "Evaluation Metrics for Labeling Comparision";
    }

    /**
     * @param labelings list of labelings to compare (restricted to 2?)
     */
    @FeatureTargetListener
    public void onLabelingsUpdated(Labeling<L>[] labelings) {
        m_labelingComparison =
                new LabelingComparison2D((NativeImgLabeling)labelings[0],
                        (NativeImgLabeling)labelings[1]);
    }

    public Img<ByteType> getDiffImg() {
        return m_labelingComparison.getDiffImg();
    }

	@Override
	public void cleanUp() {
		m_labelingComparison = null;
	}

}
