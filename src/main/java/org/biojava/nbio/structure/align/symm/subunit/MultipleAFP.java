package org.biojava.nbio.structure.align.symm.subunit;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;

public class MultipleAFP {

	AFPChain afpChain;
	Atom[] ca1;
	
	//List to store the residues aligned in the blocks. Dimensions are: [order][block_number][length of the block]
	List<List<List<Integer>>> blocks;
	//List to store the residues that are part of the free pool. Dimensions are: [order][pool_number][length of the pool]
	List<List<Integer>> free_pool;
	
	//This class pretends to use the CEMC approach for multiple structural alignment using the pairwise alignments obtained
	//from the align method.
	//TODO
	
}
