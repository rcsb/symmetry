package org.biojava.nbio.structure.align.symm;

/**
 * The internal symmetry detection can be divided into two types: 
 * CLOSED: includes the circular and dihedral symmetries, and
 * NON_CLOSED: includes the helical and protein repeats symmetries.
 * All internal symmetry cases share one property: all the subunits have the same 3D transformation.
 * 
 * A possible automatic way to discrimine between the two types is that the closed symmetry generates
 * CeSymm alignments with circular permutations (2 blocks in AFPChain), whereas the non-closed symmetry
 * generates alignments without a CP (only one block in AFPChain).
 * 
 * @author Aleix Lafita
 *
 */
public enum SymmetryType {
	CLOSED,
	NON_CLOSED;
}
