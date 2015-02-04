package org.biojava.nbio.structure.align.symm;

import org.biojava.nbio.structure.align.StructureAlignment;

import org.biojava.nbio.structure.align.ce.CeUserArgumentProcessor;

public class CeSymmUserArgumentProcessor extends CeUserArgumentProcessor{
	
	public StructureAlignment getAlgorithm() {
		return new CeSymm();
	}
}
