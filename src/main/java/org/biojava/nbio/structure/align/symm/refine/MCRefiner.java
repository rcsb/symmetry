package org.biojava.nbio.structure.align.symm.refine;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Stack;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.gui.SymmetryDisplay;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.symm.subunit.MultipleAFP;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Creates a refined alignment by a MC optimization of the subunit multiple alignment.
 * @author lafita
 */

public class MCRefiner implements Refiner {

	public MCRefiner() {
		super();
	}
	
	@Override
	public AFPChain refine(AFPChain[] afpAlignments, Atom[] ca1, Atom[] ca2, int order)
			throws RefinerFailedException {
		
		AFPChain refinedAFP = afpAlignments[0];
		
		return refinedAFP;
	}
}