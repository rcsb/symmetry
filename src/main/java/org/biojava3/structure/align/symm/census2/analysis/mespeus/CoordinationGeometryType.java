package org.biojava3.structure.align.symm.census2.analysis.mespeus;

/**
 * An ideal coordination geometry, such as trigonal planar or pentagonal bipyramidal. Also stores a type called "none" with a coordination number of 1.
 * @author dmyersturnbull
 */
public enum CoordinationGeometryType {

	NONE(1), LINEAR(2), TRIGONAL_PLANAR(3), TETRAHEDRAL(4), SQUARE_PLANAR(4), TRIGONAL_BIPYRAMIDAL(5), SQUARE_PYRAMIDAL(5),
	OCTAHEDRAL(6), TRIGONAL_PRISMATIC(6), PENTAGONAL_BIPYRAMIDAL(7), SQUARE_ANTIPRISMATIC(8), TRI_CAPPED_TRIGONAL_PRISMATIC(9);
	
	private final int coordinationNumber;

	private CoordinationGeometryType(int coordinationNumber) {
		this.coordinationNumber = coordinationNumber;
	}

	public int getCoordinationNumber() {
		return coordinationNumber;
	}
	
	
}
