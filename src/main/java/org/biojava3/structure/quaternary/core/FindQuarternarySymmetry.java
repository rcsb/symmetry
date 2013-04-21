package org.biojava3.structure.quaternary.core;

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
	private int maxFolds = 1;
	private boolean pseudoSymmetric = false;

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

	public int getChainCount() {
		run();
		return chainIds.size();
	}

	public List<String> getChainIds() {
		run();
		return chainIds;
	}

	public int getMaxFolds() {
		run();
		return maxFolds;
	}

	public Subunits getSubunits() {
		run();
		return subunits;
	}

	public String getMethod() {
		run();
		return method;
	}

	public boolean isPseudoSymmetric() {
		run();
		return pseudoSymmetric;
	}

	private void run() {
		if (modified) {
			try {
				createSubunits();
				determineRotationGroup();
				modified = false;
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}

	private void createSubunits() {
		ChainClusterer chainClusterer = new ChainClusterer(structure, parameters);
		// TODO how about chains with UNK residues??
		pseudoSymmetric = chainClusterer.isPseudoSymmetric();
		compositionFormula = chainClusterer.getCompositionFormula();
		chainIds = chainClusterer.getChainIdsInClusterOrder();
		if (chainIds.size() == 0) {
			System.err.println("createSubunits Could not find chainIds!"); 
			maxFolds = 0;
			return;
		}
		List<Integer> folds = chainClusterer.getFolds();
		maxFolds = folds.get(folds.size()-1);
		subunits = new Subunits(chainClusterer.getCalphaCoordinates(), 
				chainClusterer.getSequenceClusterIds(),
				folds,
				chainClusterer.getChainIdsInClusterOrder(),
				chainClusterer.getModelNumbersInClusterOrder());
		maxFolds = folds.get(folds.size()-1);
		if ( parameters.isVerbose())
			System.out.println("Subunits: centroids: " + subunits.getCentroid() + " folds:" + subunits.getFolds());
		//		subunits = new Subunits(chainClusterer.getCalphaCoordinates(), 
		//				chainClusterer.getSequenceClusterIds(),
		//				folds,
		//				chainClusterer.getChainIds(),
		//				chainClusterer.getModelNumbers());
		//		subunits = new Subunits(chainClusterer.getCalphaCoordinates(), 
		//				chainClusterer.getCbetaCoordinates(), 
		//				chainClusterer.getSequenceClusterIds(),
		//				folds);
	}

	private void determineRotationGroup() {
		if (chainIds.size() == 0) {
			if ( parameters.isVerbose()) {
				System.err.println("FindQuaternarySymmetry: determineRotationGroups found chainIds.size() = 0" );
			}
			return;
		}
		
		if ( parameters.isVerbose()) {
			System.err.println("FindQuaternarySymmetry: determineRotationGroups found chainIds.size() = "+ chainIds.size() );
			
			System.err.println("FindQuaternarySymmetry: determineRotationGroups getFolds.size() : "+ subunits.getFolds().size());
			System.err.println("FindQuaternarySymmetry: determineRotationGroups subunit count: "+ subunits.getSubunitCount());
			
		}
		
		
		
		if (subunits.getFolds().size() == 1) {
			
			// no symmetry possible, create empty ("C1") rotation group
			method = "norotation";
			
			if ( parameters.isVerbose()) {
				System.out.println("FindQuaternarySymmetry: determineRotationGroups: " + method);
			}
			
			rotationGroup =  new RotationGroup();
			rotationGroup.setC1(subunits.getSubunitCount());
		} else if (subunits.getSubunitCount() == 2 && subunits.getFolds().contains(2)) {
			method = "C2rotation";
			
			if ( parameters.isVerbose()) {
				System.out.println("FindQuaternarySymmetry: determineRotationGroups: " + method);
			}
			
			QuatSymmetrySolver solver = new C2RotationSolver(subunits, parameters);
			rotationGroup = solver.getSymmetryOperations();
		} else {
			method = "rotation";
		
			if ( parameters.isVerbose()) {
				System.out.println("FindQuaternarySymmetry: determineRotationGroups: " + method);
			}
			
			QuatSymmetrySolver solver = new RotationSolver(subunits, parameters);
			// TODO
			//			QuatSymmetrySolver solver = new SystematicSolver(subunits, parameters);
			rotationGroup = solver.getSymmetryOperations();
		}
		// should check in w/ .complete() to make sure this is called under all circumstances
		rotationGroup.complete();
	//	HelixCheck hc = new HelixCheck(subunits, rotationGroup, this.parameters);
	//	System.out.println("Helical: " + hc.isHelical());
	//	System.exit(-1);
	}
}
