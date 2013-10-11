package org.biojava3.structure.align.symm.census2;

import org.biojava.bio.structure.align.model.AFPChain;

/**
 * An alignment function.
 * In future versions, may be able to handle different mappings between different symmetry subunits.
 * Also critical for serializing AFPChains within Result XML files.
 * @author dmyersturnbull
 */
public class AlignmentMapping {

	private AFPChain afpChain;
	
	public AlignmentMapping(AFPChain afpChain) {
		this.afpChain = afpChain;
	}

	public AFPChain toAfpChain() {
		return afpChain;
	}

}
