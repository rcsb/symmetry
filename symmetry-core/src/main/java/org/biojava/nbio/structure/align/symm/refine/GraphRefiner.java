package org.biojava.nbio.structure.align.symm.refine;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.symmetry.internal.Refiner;
import org.biojava.nbio.structure.symmetry.internal.RefinerFailedException;

/**
 * This refinement transforms the self-alignment into a Graph and
 * extracts its Components. It then generates a MultipleAlignment
 * by comibining the connected Components.
 * 
 * @author Aleix Lafita
 * 
 */
public class GraphRefiner implements Refiner {

	@Override
	public AFPChain refine(List<AFPChain> afpAlignments, Atom[] atoms,
			int order) throws StructureException, RefinerFailedException {
		return null;
	}
}
