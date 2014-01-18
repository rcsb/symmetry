package org.biojava3.structure.align.symm.order;

import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;

/**
 * Calls Spencer's method for determining order.
 * @author dmyersturnbull
 */
public class SequenceFunctionOrderDetector implements OrderDetector {

	private int maxSymmetry = 8;
	private float minimumMetricChange = 0.4f;
	
	public SequenceFunctionOrderDetector() {
		super();
	}

	public SequenceFunctionOrderDetector(int maxSymmetry, float minimumMetricChange) {
		super();
		this.maxSymmetry = maxSymmetry;
		this.minimumMetricChange = minimumMetricChange;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {
		try {
			Map<Integer,Integer> alignment = AlignmentTools.alignmentAsMap(afpChain);

			return AlignmentTools.getSymmetryOrder(alignment,
					new AlignmentTools.IdentityMap<Integer>(), maxSymmetry, minimumMetricChange);
//			return AlignmentTools.getSymmetryOrder(afpChain, maxSymmetry, minimumMetricChange);
		} catch (StructureException e) {
			throw new OrderDetectionFailedException(e);
		}
	}

}
