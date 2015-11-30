package org.biojava.nbio.structure.align.symm.order;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.OrderDetector;
import org.biojava.nbio.structure.symmetry.internal.RefinerFailedException;

public class MultipassOrderDetector implements OrderDetector {

	private int maxOrder = 8;

	public MultipassOrderDetector() {
		this(8);
	}
	/**
	 */
	public MultipassOrderDetector(int maxOrder) {
		this.maxOrder = maxOrder;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca)
			throws RefinerFailedException {
		CeSymm ce = new CeSymm();
		CESymmParameters params = (CESymmParameters) ce.getParameters();
		params.setMaxSymmOrder(maxOrder);
		params.setRefineMethod(RefineMethod.MULTIPLE);
		params.setSymmType(SymmetryType.CLOSE);
		params.setOptimization(false);
		try {
			ce.analyze(ca);
		} catch (StructureException e) {
			throw new RefinerFailedException(e);
		}
		List<AFPChain> alignments = ce.getSelfAlignments();
		// For high orders, take it from the number of alignments with unrefined TM above threshold
		if(alignments.size() > 1) {
			return alignments.size() + 1;
		}
		// For C1/C2, look at refined TM
		if(alignments.get(0).getTMScore() >= 
				CESymmParameters.DEFAULT_SYMMETRY_THRESHOLD){
			return 2;
		} else {
			return 1;
		}
	}
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "MultipassOrderDetector [maxOrder=" + maxOrder + "]";
	}

}
