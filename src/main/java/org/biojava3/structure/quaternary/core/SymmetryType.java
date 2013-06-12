package org.biojava3.structure.quaternary.core;

import java.io.Serializable;

public enum SymmetryType implements Serializable{
	Asymmetric("Asymmetric"), 
	GlobalSymmetry("GlobalSymmetry"), 
	GlobalPseudoSymmetry("GlobalPseudoSymmetry"),
	LocalSymmetry("LocalSymmetry"), 
	GlobalPseudosymmetrStructureOnly("GlobalPseudosymmetrStructureOnly"), 
	LocalPseudosymmetryStructureOnly("LocalPseudosymmetryStructureOnly");

	public final String type;
	
	SymmetryType(String type){
		this.type = type;
	}
	
	public static SymmetryType getSymmetryTypeFromString(String type){
		 for(SymmetryType rt : SymmetryType.values()) {
			 if (rt.type.equalsIgnoreCase(type)) {
				 return rt;
			 }
		 }
		 
		 return null;
	}
}