package org.biojava.nbio.structure.align.symm.refine;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * A method to refine the AFP alignment from one or more self-alignments in order to make the subunits consistent.
 * @author lafita
 *
 */
public interface Refiner {

	AFPChain refine(AFPChain[] afpAlignments, Atom[] ca1, Atom[] ca2, int order) throws RefinerFailedException;
	
}
