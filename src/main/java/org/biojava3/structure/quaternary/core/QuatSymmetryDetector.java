/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2013-05-23
 *
 */
package org.biojava3.structure.quaternary.core;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.Structure;
import org.biojava3.structure.quaternary.misc.CombinationGenerator;
import org.biojava3.structure.quaternary.utils.ComponentFinder;
import org.biojava3.structure.quaternary.utils.Graph;

/**
 * Detects global and local quaternary protein structure symmetry in a structure.
 * 
 * The QuatSymmetryParameter settings affect the calculated results.
 * 
 * @author Peter Rose
 *
 */
public class QuatSymmetryDetector {
	private Structure structure = null;
	private QuatSymmetryParameters parameters = null;
	
	private ChainClusterer chainClusterer = null;
	private List<QuatSymmetryResults> globalSymmetry = new ArrayList<QuatSymmetryResults>();
	private List<List<QuatSymmetryResults>> localSymmetries = new ArrayList<List<QuatSymmetryResults>>();
	private boolean complete = false;

	public QuatSymmetryDetector(Structure structure, QuatSymmetryParameters parameters) {
		this.structure = structure;
		this.parameters = parameters;
	}
	
	/**
	 * Returns true if structure contains protein subunits. The other methods
	 * in this class will only return data if protein subunits are present.
	 * Always use this method first, before retrieving global or local symmetry 
	 * results.
	 * 
	 * @return true if protein subunits are present
	 */
	public boolean hasProteinSubunits() {
		run();
		return chainClusterer.getChainIds().size() > 0;
	}
	
	/**
	 * Returns list of global quaternary structure symmetry results
	 * 
	 * @return list of global quaternary structure symmetry results
	 */
	public List<QuatSymmetryResults> getGlobalSymmetry() {
		run();
		return globalSymmetry;
	}
	
	/**
	 * Returns a list of lists of local quaternary structure symmetry results
	 * 
	 * @return list of lists of local quaternary structure symmetry results
	 */
	public List<List<QuatSymmetryResults>> getLocalSymmetries() {
		run();
		return localSymmetries;
	}
	
	private void run() {
		if (complete) {
			return;
		}
		complete = true;
		ClusterProteinChains clusterer = new ClusterProteinChains(structure, parameters);
		
		// sort seq. identity thresholds from smallest to largest. This reduces the total number of calculations necessary.
		double[] thresholds = parameters.getSequenceIdentityThresholds().clone();
		Arrays.sort(thresholds);
		
		for (int index = 0; index < thresholds.length; index++) {
			chainClusterer = new ChainClusterer(clusterer.getSequenceAlignmentClusters(thresholds[index]));

			int chainCount = chainClusterer.getChainIds().size();
			int clusterCount = chainClusterer.getSequenceClusterCount();

			if (chainCount == 0) {
				return;
			}

			// determine global symmetry
			Subunits globalSubunits = createGlobalSubunits();
			QuatSymmetryResults gSymmetry = calcQuatSymmetry(globalSubunits);
			gSymmetry.setSequenceIdentityThreshold(thresholds[index]);
			globalSymmetry.add(gSymmetry);

			// determine local symmetry if global structure is 
			// (1) asymmetric (C1)
			// (2) heteromeric (belongs to more than 1 sequence cluster)
			// (3) more than 2 chains (heteromers with just 2 chains cannot have local symmetry)
			
			// TODO example 2PT7: global C2, but local C6 symm., should that be included here ...?
			// i.e., include all heteromers here, for example if higher symmetry is possible by stoichiometry? A6B2 -> local A6  can have higher symmetry
			if (parameters.isLocalSymmetry()) {
				if (gSymmetry.getRotationGroup().getPointGroup().equals("C1") &&
						clusterCount > 1 && chainCount > 2) {

					List<Subunits> localSubunits = createLocalSubunits();
					List<QuatSymmetryResults> lSymmetry = new ArrayList<QuatSymmetryResults>();

					for (Subunits subunits: localSubunits) {
						QuatSymmetryResults result = calcQuatSymmetry(subunits);
						addToLocalSymmetry(result, lSymmetry);
					}
					localSymmetries.add(lSymmetry);
				}
			}
			
			if (! gSymmetry.getSubunits().isPseudoStoichiometric()) {
				break;
			}
		}
		
		setPreferredResults();
		setPseudoSymmetry();
	}
	
	/**
	 * Sets preferred results flag for symmetry result that should be shown by default in visualization programs
	 * @param thresholds sequence identity thresholds
	 */
	private void setPreferredResults() {
		int[] score = new int[globalSymmetry.size()];
		
		// check global symmetry
		for (int i = 0; i < globalSymmetry.size(); i++) {
			QuatSymmetryResults result = globalSymmetry.get(i);
			if (! result.getRotationGroup().getPointGroup().equals("C1")) {
				score[i] += 2;
			}
			if (! result.getSubunits().isPseudoStoichiometric()) {
				score[i]++;
			}
		}

		int bestGlobal = 0;
		int bestScore = 0;
		for (int i = 0; i < score.length; i++) {
			if (score[i] > bestScore) {
				bestScore = score[i];
				bestGlobal = i;
			}
		}
		if (bestScore >= 2) {
			QuatSymmetryResults g = globalSymmetry.get(bestGlobal);
			g.setPreferredResult(true);
			return;
		}

		// check local symmetry
		score = new int[localSymmetries.size()];

		for (int i = 0; i < localSymmetries.size(); i++) {
			List<QuatSymmetryResults> results = localSymmetries.get(i);
			for (QuatSymmetryResults result: results) {
				if (! result.getRotationGroup().getPointGroup().equals("C1")) {
					score[i] += 2;
				}
				if (! result.getSubunits().isPseudoStoichiometric()) {
					score[i]++;
				}
			}
		}
	
		int bestLocal = 0;
		bestScore = 0;
		for (int i = 0; i < score.length; i++) {
			if (score[i] > bestScore) {
				bestScore = score[i];
				bestLocal = i;
			}
		}
		if (bestScore > 0) {
			List<QuatSymmetryResults> results = localSymmetries.get(bestLocal);
			for (QuatSymmetryResults result: results) {
				result.setPreferredResult(true);
			}
		} else {
			QuatSymmetryResults g = globalSymmetry.get(bestGlobal);
			g.setPreferredResult(true);
		}
	}
	
	/**
	 * Sets pseudosymmetry flag for results that have pseudosymmetry
	 * @param thresholds sequence identity thresholds
	 */
	private void setPseudoSymmetry() {
		String symmPointGroup = "";
		String pseudoPointGroup = "";
		QuatSymmetryResults pseudoGlobal = null;
		for (QuatSymmetryResults result: getGlobalSymmetry()) {
			if (result.getSubunits().isPseudoStoichiometric()) {
				pseudoPointGroup = result.getRotationGroup().getPointGroup();
				pseudoGlobal = result;
			} else {
				symmPointGroup = result.getRotationGroup().getPointGroup();
			}
		}

		if (pseudoGlobal != null && ! pseudoPointGroup.equals(symmPointGroup)) {
			pseudoGlobal.getSubunits().setPseudoSymmetric(true);
		}

		// check for local pseudosymmetry
		List<QuatSymmetryResults> pseudoLocal = null;
		for (List<QuatSymmetryResults> results: getLocalSymmetries()) {
			for (QuatSymmetryResults result: results) {
				if (result.getSubunits().isPseudoStoichiometric()) {
					pseudoPointGroup = result.getRotationGroup().getPointGroup();
					pseudoLocal = results;
				} else {
					symmPointGroup = result.getRotationGroup().getPointGroup();
				}
			}
		}
		if (pseudoLocal != null && ! pseudoPointGroup.equals(symmPointGroup)) {
			for (QuatSymmetryResults result: pseudoLocal) {
				result.getSubunits().setPseudoSymmetric(true);
			}
		}
	}

	private void addToLocalSymmetry(QuatSymmetryResults testResults, List<QuatSymmetryResults> localSymmetry) {
		if (testResults.getRotationGroup().getPointGroup().equals("C1")) {
			return;
		}

		for (QuatSymmetryResults results: localSymmetry) {
			if (results.getSubunits().overlaps(testResults.getSubunits())) {
//				System.out.println("Subunit overlaps large subunit");
				if (testResults.getRotationGroup().getOrder() <= results.getRotationGroup().getOrder()) {
//					System.out.println("Subunit has equal or lower symmetry order");
					return;
				}
			}
		}
		testResults.setLocal(true);
		localSymmetry.add(testResults);
	}

	private Subunits createGlobalSubunits() {
		Subunits subunits = new Subunits(chainClusterer.getCalphaCoordinates(), 
				chainClusterer.getSequenceClusterIds(),
				chainClusterer.getPseudoStoichiometry(),
				chainClusterer.getMinSequenceIdentity(),
				chainClusterer.getMaxSequenceIdentity(),
				chainClusterer.getFolds(),
				chainClusterer.getChainIds(),
				chainClusterer.getModelNumbers());
		return subunits;
	}
	
	private List<Subunits> createLocalSubunits() {
		List<Subunits> subunits = new ArrayList<Subunits>();
		List<List<Integer>> subClusters = decomposeClusters(chainClusterer.getCalphaCoordinates(), chainClusterer.getSequenceClusterIds());
		for (List<Integer> subCluster: subClusters) {
			subunits.add(createLocalSubunit(subCluster));
		}
		return subunits;
	}
	
	private Subunits createLocalSubunit(List<Integer> subCluster) {
	      List<Point3d[]> subCalphaCoordinates = new ArrayList<Point3d[]>(subCluster.size());   
	      List<Integer> subSequenceIds = new ArrayList<Integer>(subCluster.size());
	      List<Boolean> subPseudoStoichiometry = new ArrayList<Boolean>(subCluster.size());
	      List<Double> subMinSequenceIdentity = new ArrayList<Double>(subCluster.size());
	      List<Double> subMaxSequenceIdentity = new ArrayList<Double>(subCluster.size());
	      List<String> subChainIds = new ArrayList<String>(subCluster.size());
	      List<Integer> subModelNumbers = new ArrayList<Integer>(subCluster.size());

	      for (int index: subCluster) {
	    	  subCalphaCoordinates.add(chainClusterer.getCalphaCoordinates().get(index));
	    	  subSequenceIds.add(chainClusterer.getSequenceClusterIds().get(index));
	    	  subPseudoStoichiometry.add(chainClusterer.getPseudoStoichiometry().get(index));
	    	  subMinSequenceIdentity.add(chainClusterer.getMinSequenceIdentity().get(index));
	    	  subMaxSequenceIdentity.add(chainClusterer.getMaxSequenceIdentity().get(index));
	    	  subChainIds.add(chainClusterer.getChainIds().get(index));
	    	  subModelNumbers.add(chainClusterer.getModelNumbers().get(index));
	      }

	      standardizeSequenceIds(subSequenceIds);
	      
	      Integer[] array = subSequenceIds.toArray(new Integer[subSequenceIds.size()]);
	      List<Integer> subFolds = getFolds(array, subSequenceIds.size());
	      Subunits subunits = new Subunits(subCalphaCoordinates, 
					subSequenceIds,
					subPseudoStoichiometry,
					subMinSequenceIdentity,
					subMaxSequenceIdentity,
			        subFolds,
					subChainIds,
					subModelNumbers);
			return subunits;
	}

	/**
	 * Resets list of arbitrary sequence ids into integer order: 0, 1, ...
	 * @param subSequenceIds
	 */
	private static void standardizeSequenceIds(List<Integer> subSequenceIds) {
		int count = 0;
	      int current = subSequenceIds.get(0);
	      for (int i = 0; i < subSequenceIds.size(); i++) {
	    	  if (subSequenceIds.get(i) > current) {
	    		  current = subSequenceIds.get(i);
	    		  count++;
	    	  }
	    	  subSequenceIds.set(i, count);
	      }
	}
	
	private List<List<Integer>> decomposeClusters(List<Point3d[]> caCoords, List<Integer> clusterIds) {
		List<List<Integer>> subClusters = new ArrayList<List<Integer>>();

		int last = getLastMultiSubunit(clusterIds);
		List<Point3d[]> subList = caCoords;
		if (last < caCoords.size()) {
			subList = caCoords.subList(0, last);
		} else {
			last = caCoords.size();
		}

		SubunitGraph subunitGraph = new SubunitGraph(subList);
		Graph<Integer> graph = subunitGraph.getProteinGraph();
//		System.out.println("Graph: " + graph);

		for (int i = last; i > 1; i--) {
			CombinationGenerator generator = new CombinationGenerator(last, i);
			int[] indices = null;
			Integer[] subCluster = new Integer[i];
			
			// avoid combinatorial explosion, i.e. for 1FNT
			BigInteger maxCombinations = BigInteger.valueOf(parameters.getMaximumLocalCombinations());
		    if (generator.getTotal().compareTo(maxCombinations) > 0) {
		    	continue;
		    }
			
			while (generator.hasNext()) {
				indices = generator.getNext();	
				List<Integer> subSet = new ArrayList<Integer>(indices.length);
				for (int index: indices) {
					subSet.add(index);
				}
				Graph<Integer> subGraph = graph.extractSubGraph(subSet);		
				//			System.out.println("Subgraph: " + subGraph);

				if (isConnectedGraph(subGraph)) {		
					for (int j = 0; j < indices.length; j++) {
						subCluster[j] = clusterIds.get(indices[j]);
					}		
					List<Integer> folds = getFolds(subCluster, last);
					if (folds.size() > 1) {
						subClusters.add(subSet);
					}
				}
			}
		}

		return subClusters;
	}

	private static int getLastMultiSubunit(List<Integer> clusterIds) {
		for (int i = 0, n = clusterIds.size(); i < n; i++) {
			if (i < n-2) {
				if (clusterIds.get(i)!=clusterIds.get(i+1) && 
						clusterIds.get(i+1) != clusterIds.get(i+2)) {
					return i+1;
				}
			}
			if (i == n-2) {
				if (clusterIds.get(i)!=clusterIds.get(i+1)) {
					return i+1;
				}
			}
		}
		return clusterIds.size();
	}
	
	private static boolean isConnectedGraph(Graph<Integer> graph) {
		ComponentFinder<Integer> finder = new ComponentFinder<Integer>();
		finder.setGraph(graph);
		return finder.getComponentCount() == 1;
	}
	
	private static List<Integer> getFolds(Integer[] subCluster, int size) {
		List<Integer> denominators = new ArrayList<Integer>();
		int[] counts = new int[size];
		for (int element: subCluster) {
			counts[element]++;
		}

		for (int d = 1; d <= subCluster.length; d++) {
			int count = 0;
			for (int i = 0; i < size; i++) {
				if (counts[i] > 0 && (counts[i] % d == 0)) {
					count += counts[i];
				}
			}
			if (count == subCluster.length) {
				denominators.add(d);
			}
		}
		
		Collections.sort(denominators);
		return denominators;
	}
	
	private QuatSymmetryResults calcQuatSymmetry(Subunits subunits) {
		if (subunits.getSubunitCount() == 0) {
			return null;
		}
		
		RotationGroup rotationGroup = null;
		String method = null;
		if (subunits.getFolds().size() == 1) {			
			// no symmetry possible, create empty ("C1") rotation group
			method = "norotation";
			rotationGroup =  new RotationGroup();
			rotationGroup.setC1(subunits.getSubunitCount());
		} else if (subunits.getSubunitCount() == 2 && subunits.getFolds().contains(2)) {
			method = "C2rotation";
			QuatSymmetrySolver solver = new C2RotationSolver(subunits, parameters);
			rotationGroup = solver.getSymmetryOperations();
		} else {
			method = "rotation";
			QuatSymmetrySolver solver = new RotationSolver(subunits, parameters);
			rotationGroup = solver.getSymmetryOperations();
		}

		rotationGroup.complete();
		
//		String pointGroup = rotationGroup.getPointGroup();
//		if (pointGroup.startsWith("C")) {
//			HelixCheck hc = new HelixCheck(subunits, rotationGroup, this.parameters);
//			System.out.println("Helical: " + hc.isHelical());
//		}
	//	System.exit(-1);
		QuatSymmetryResults results = new QuatSymmetryResults(subunits, rotationGroup, method);
		return results;
	}
}
