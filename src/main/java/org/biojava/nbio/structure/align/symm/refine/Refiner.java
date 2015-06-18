package org.biojava.nbio.structure.align.symm.refine;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * A method to refine the AFP alignment from one or more alternative self-alignments in order to make the subunits consistent.
 * @author Aleix Lafita
 *
 */
public interface Refiner {

	/**
	 * Returns a refined symmetry alignment, where the subunit residues are aligned consistently and separated
	 * into the blocks of the AFPChain.
	 * 
	 * @param afpAlignments List of returned self-alignments in CeSymm
	 * @param ca1 coordinates of the structure
	 * @param ca2 coordinates of the structure
	 * @param order of symmetry detected by the OrderDetector
	 * @return AFPChain refined symmetry alignment
	 * @throws RefinerFailedException
	 * @throws StructureException
	 */
	AFPChain refine(List<AFPChain> afpAlignments, Atom[] ca1, Atom[] ca2, int order) throws RefinerFailedException,StructureException;
	
}
