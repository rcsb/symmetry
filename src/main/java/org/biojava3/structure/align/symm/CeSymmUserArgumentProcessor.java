package org.biojava3.structure.align.symm;

import org.biojava.bio.structure.align.StructureAlignment;

import org.biojava.bio.structure.align.ce.CeUserArgumentProcessor;

public class CeSymmUserArgumentProcessor extends CeUserArgumentProcessor{
	
	public StructureAlignment getAlgorithm() {
		return new CeSymm();
	}
}
