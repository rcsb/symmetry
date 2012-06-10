package org.biojava3.structure.align.symm.quaternary;

import java.util.List;

import org.biojava.bio.structure.Structure;

public class FindQuarternarySymmetry {
	private Structure structure = null;
	private QuatSymmetryParameters parameters = null;

	private RotationGroup rotationGroup = null;
	private String compositionFormula = "";
	private List<String> chainIds = null;
	private Subunits subunits = null;
	private String method = "";
	
	private boolean modified = true;

	public FindQuarternarySymmetry(Structure structure, QuatSymmetryParameters parameters) {
		this.structure = structure;
		this.parameters = parameters;
	}

	public RotationGroup getRotationGroup() {
		run();
        return rotationGroup;
	}
	
	public String getCompositionFormula() {
		run();
		return compositionFormula;
	}
	
	public List<String> getChainIds() {
		run();
		return chainIds;
	}
	
	public Subunits getSubunits() {
		run();
		return subunits;
	}
	
	public String getMethod() {
		run();
		return method;
	}
	
	private void run() {
		if (modified) {
			createSubunits();
			determineRotationGroup();
			modified = false;
		}
	}

	private void createSubunits() {
		ChainClusterer chainClusterer = new ChainClusterer(structure, parameters);
		// TODO how about chains with UNK residues??
		compositionFormula = chainClusterer.getCompositionFormula();
		chainIds = chainClusterer.getChainIds();
		subunits = new Subunits(chainClusterer.getCalphaCoordinates(), chainClusterer.getCbetaCoordinates(), chainClusterer.getSequenceClusterIds());
	}

	private void determineRotationGroup() {
		if (subunits.getSubunitCount() <= 1) {
			rotationGroup =  new RotationGroup();
			method = "none";
			return;
		} else if (subunits.getSubunitCount() == 2) {
			method = "C2rotation";
			QuatSymmetrySolver solver = new C2RotationSolver(subunits, parameters.getRmsdThreshold());
			rotationGroup = solver.getSymmetryOperations();
		} else {
			System.out.println("Rotation solver");
			method = "rotation";
			QuatSymmetrySolver solver = new RotationSolver(subunits, parameters.getRmsdThreshold());
			rotationGroup = solver.getSymmetryOperations();
		}
		// TODO use composition formula to calculate the maximum number of possible symmetry operations and check here -> try sys. solver?
	}
}
