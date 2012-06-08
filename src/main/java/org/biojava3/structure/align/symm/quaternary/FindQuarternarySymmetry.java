package org.biojava3.structure.align.symm.quaternary;

import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.Structure;

public class FindQuarternarySymmetry {
	private Structure structure = null;
	private List<Point3d[]> caCoords = null;
	private List<Point3d[]> cbCoords = null;
	private List<Integer> sequenceClusterIds = null;
	private boolean pseudoSymmetryAllowed = false;
	private Subunits subunits = null;
	private ChainClusterer chainClusterer = null;
	private RotationGroup symmetryOperations = null;
	private String method = "";
	private int minimumSequenceLength = 24;
	private double sequenceIdentityThreshold = 1.0;

	public FindQuarternarySymmetry(Structure structure) {
		this.structure = structure;
	}
	
	public void setPseudoSymmetryAllowed(boolean pseudoSymmetryAllowed) {
		this.pseudoSymmetryAllowed = pseudoSymmetryAllowed;
	}

	public RotationGroup getRotationGroup() {
		run();
        return symmetryOperations;
	}
	
	public Subunits getSubunits() {
		return subunits;
	}
	
	public ChainClusterer getSequenceCluster() {
		return chainClusterer;
	}
	
	public int getChainCount() {
		return caCoords.size();
	}
	
	public String getMethod() {
		return method;
	}
	
	public void setMinimumSequenceLength(int minimumSequenceLength) {
		this.minimumSequenceLength = minimumSequenceLength;
	}
	
	private boolean run() {
		chainClusterer = new ChainClusterer(structure, minimumSequenceLength, sequenceIdentityThreshold);
		String formula = chainClusterer.getCompositionFormula();
		System.out.println("Formula: " + formula);
		System.out.println(chainClusterer);

		caCoords = chainClusterer.getCalphaCoordinates();
		cbCoords = chainClusterer.getCbetaCoordinates();
		sequenceClusterIds = chainClusterer.getSequenceClusterIds();
	    subunits = new Subunits(caCoords, cbCoords, sequenceClusterIds);
	    QuatSymmetryPerceptor perceptor = new QuatSymmetryPerceptor(subunits);
	    perceptor.setPseudoSymmetryAllowed(pseudoSymmetryAllowed);
	    symmetryOperations = perceptor.getSymmetryOperations();
	    method = perceptor.getMethod();

	           System.out.println("--- SymmetryOperations ---");;
	           System.out.println(symmetryOperations);
	    	              System.out.println(symmetryOperations.getPointGroup());
	    return true;
	}
}
