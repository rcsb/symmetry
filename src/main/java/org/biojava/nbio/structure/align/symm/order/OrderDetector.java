package org.biojava.nbio.structure.align.symm.order;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * A method to decide the order of symmetry given a self-alignment.
 * @author dmyersturnbull
 *
 */
public interface OrderDetector {

	int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException;
	
}
