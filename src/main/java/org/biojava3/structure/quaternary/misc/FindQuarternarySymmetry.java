package org.biojava3.structure.quaternary.misc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.Structure;
import org.biojava3.structure.quaternary.core.C2RotationSolver;
import org.biojava3.structure.quaternary.core.ChainClusterer;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.QuatSymmetrySolver;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.core.RotationSolver;
import org.biojava3.structure.quaternary.core.SubunitGraph;
import org.biojava3.structure.quaternary.core.Subunits;
import org.biojava3.structure.quaternary.utils.Graph;

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
		compositionFormula = chainClusterer.getStoichiometry();
		chainIds = chainClusterer.getChainIds();

		if (chainIds.size() == 0) {
			System.err.println("createSubunits Could not find chainIds!"); 
			maxFolds = 0;
			return;
		}
		List<Integer> folds = chainClusterer.getFolds();
		maxFolds = folds.get(folds.size()-1);
		subunits = new Subunits(chainClusterer.getCalphaCoordinates(), 
				chainClusterer.getSequenceClusterIds(),
				chainClusterer.getPseudoStoichiometry(),
				folds,
				chainClusterer.getChainIds(),
				chainClusterer.getModelNumbers());
		pseudoSymmetric = subunits.isPseudoStiochiometric();
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
	
	private void decomposeSubunits(List<Point3d[]> caCoords, List<Integer> clusterIds) {
		System.out.println("ClusterIds: " + clusterIds);
		SubunitGraph subunitGraph = new SubunitGraph(caCoords);
		Graph<Integer> graph = subunitGraph.getProteinGraph();
		for (int i = 2; i < graph.size()+1; i++) {
			CombinationGenerator generator = new CombinationGenerator(graph.size(), i);
			int[] indices = new int[i];
			int[] subCluster = new int[i];
			while (generator.hasNext()) {
				indices = generator.getNext();
				
				if (isConnectedGraph(graph, indices)) {
					for (int j = 0; j < indices.length; j++) {
						subCluster[j] = clusterIds.get(indices[j]);
					}		
					List<Integer> folds = getFolds(subCluster, graph.size());
					if (! folds.isEmpty()) {
						System.out.println("Complex:    " + Arrays.toString(indices));
						System.out.println("Subcluster: " + Arrays.toString(subCluster));
						System.out.println("Folds:      " + folds);
					}
				}
			}
		}
		System.out.println("Subunit graph");
		System.out.println(graph);
	}
	
	private boolean isConnectedGraph(Graph<Integer> graph, int[] indices) {
		for (int i = 0; i < indices.length; i++) {
			boolean connected = false;
			for (int j = 1; j < indices.length; j++) {
				if (i == j) continue;
				if (graph.containsEdge(indices[i], indices[j])) {
					connected = true;
					break;
				}
			}
			if (!connected) {
				return false;
			}
		}
		return true;
	}
	
	private List<Integer> getFolds(int[] clusters, int size) {
		List<Integer> denominators = new ArrayList<Integer>();
		int[] counts = new int[size];
		for (int element: clusters) {
			counts[element]++;
		}
//		System.out.println("Counts: " + Arrays.toString(counts));
		for (int d = 2; d <= clusters.length; d++) {
			int count = 0;
			for (int i = 0; i < size; i++) {
				if (counts[i] > 0 && (counts[i] % d == 0)) {
	//				count++;
					count += counts[i];
				}
			}
//			System.out.println("d, count: " + d + ": " + count);
			if (count == clusters.length) {
				denominators.add(d);
			}
		}
		
		Collections.sort(denominators);
		return denominators;
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
		
//		String pointGroup = rotationGroup.getPointGroup();
//		if (pointGroup.startsWith("C")) {
//			HelixCheck hc = new HelixCheck(subunits, rotationGroup, this.parameters);
//			System.out.println("Helical: " + hc.isHelical());
//		}
	//	System.exit(-1);
	}
}
