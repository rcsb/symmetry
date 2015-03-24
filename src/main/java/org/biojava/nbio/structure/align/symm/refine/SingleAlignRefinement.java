package org.biojava.nbio.structure.align.symm.refine;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.SymmRefiner;

/**
 * Creates a refined alignment with a single self-alignment. Needs the order of symmetry.
 * Uses Spencer's refinement method.
 * @author lafita
 */

public class SingleAlignRefinement implements Refiner {

	public SingleAlignRefinement() {
		super();
	}
	
	@Override
	public AFPChain refine(AFPChain[] afpAlignments, Atom[] ca1, Atom[] ca2, int order)
			throws RefinerFailedException {
		
		AFPChain originalAFP = afpAlignments[0];
		AFPChain refinedAFP = new AFPChain();
		
		//Provisional, move the code from SymRefiner to this class.
		try {
			refinedAFP = SymmRefiner.refineSymmetry(originalAFP, ca1, ca2, order);
		} catch (StructureException e) {
			e.printStackTrace();
		}
		
		return refinedAFP;
	}

}
