package org.biojava.nbio.structure.align.symm;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.gui.SymmetryDisplay;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class provides methods for sorting the chains  of a structure 
 * in their symmetric order, so that the CeSymm algorithm (and other 
 * sequence-dependent algorithms) can deal with multiple chains.
 * 
 * @author Aleix Lafita
 *
 */
public class ChainSorter {

	private static final Logger logger = 
			LoggerFactory.getLogger(ChainSorter.class);
	
	/**
	 * Application: Cyclic Symmetry (Cn).
	 * <p>
	 * This method sorts the chains by the following steps:
	 * <ul><li>Calculate the centroid of each chain.
	 * <li>Pick the farthest chain from all others as the first 
	 * 		(for open cases, can be random for closed).
	 * <li>Iteratively pick the closest chain to the last one from 
	 * 		the remaining set.
	 * </ul>
	 * In case of homo n-meric cyclic symmetry (Cn) the following 
	 * theorem applies: "the closest point of another point is either 
	 * the previous or the next in the 3D rotation around the axis of 
	 * symmetry". That is true because the points are situated at
	 * the edge of a circle. Thus, this property can be used to sort 
	 * the chains in cyclic order.
	 * 
	 * @param structure Structure containing the Chains. If multiple
	 * 			models are present and the structure is not NMR, the
	 * 			Chains in all models will be considered.
	 * @return Atom[] with the sorted order of Atoms, 
	 * 			corresponding to the ordering of the chains.
	 */
	public static Atom[] cyclicSort(Structure structure){
		
		List<Atom[]> chainAtoms = new ArrayList<Atom[]>();
		//Consider all models of the structure
		for (int m=0; m < structure.nrModels(); m++){
			for (Chain c:structure.getChains(m)){
				Atom[] atoms = StructureTools.getRepresentativeAtomArray(c);
				if (atoms.length > 0) chainAtoms.add(atoms);
			}
			if (structure.isNmr()) break; //Only consider first model
		}
		return cyclicSort(chainAtoms);
	}
	
	/**
	 * Application: Cyclic Symmetry (Cn).
	 * <p>
	 * This method sorts the chains by the following steps:
	 * <ul><li>Calculate the centroid of each chain.
	 * <li>Pick the farthest chain from all others as the first 
	 * 		(for open cases, can be random for closed).
	 * <li>Iteratively pick the closest chain to the last one from 
	 * 		the remaining set.
	 * </ul>
	 * In case of homo n-meric cyclic symmetry (Cn) the following 
	 * theorem applies: "the closest point of another point is either 
	 * the previous or the next in the 3D rotation around the axis of 
	 * symmetry". That is true because the points are situated at
	 * the edge of a circle. Thus, this property can be used to sort 
	 * the chains in cyclic order.
	 * 
	 * @param chains List of Atoms of each chain
	 * @return Atom[] with the sorted order of Atoms, corresponding 
	 * 			to the ordering of the chains.
	 */
	public static Atom[] cyclicSort(List<Atom[]> chains){
		
		logger.info("Ordering of "+chains.size()+" chains:");
		//Initialize the variables nedeed
		List<Point3d> centroids = new ArrayList<Point3d>();
		Matrix chainDist = new Matrix(chains.size(), chains.size());
		List<Atom> sortedAtoms = new ArrayList<Atom>();
		List<Integer> remainingChains = new ArrayList<Integer>();
		for (int n=0; n<chains.size(); n++) 
			remainingChains.add(n); //Contains the indices of the unseen chains
		
		/** STEP 1: Centroid calculation */
		//Calculate the centroids of the chains
		for (Atom[] array:chains){
			Atom centroid = Calc.getCentroid(array);
			centroids.add(new Point3d(centroid.getCoords()));
		}
		
		/** STEP 2: Pick the farthest chain to all others */
		//Calculate the interchain chain distances
		for (int m=0; m<centroids.size(); m++){
			for (int n=0; n<centroids.size(); n++){
				chainDist.set(m, n, 
						centroids.get(m).distance(centroids.get(n)));
			}
		}
		//Pick the farthest point to all others
		double maxDist = 0.0;
		int maxIndex = 0;
		for (int p=0; p<centroids.size(); p++){
			double distance = 0.0;
			for (int p2=0; p2<centroids.size(); p2++) {
				distance += chainDist.get(p, p2);
			}
			if (distance>maxDist){
				maxIndex = p;
				maxDist = distance;
			}
		}
		
		/** STEP 3: Iteratively pick the closest chain to the last one */
		//First add the farthest chain to the beginning of the array
		for (Atom a:chains.get(maxIndex)) sortedAtoms.add(a);
		remainingChains.remove((Integer) maxIndex);
		logger.info(chains.get(maxIndex)[0].getGroup().getChainId());
		
		int lastIndex = maxIndex;
		while (remainingChains.size()>0) {
			double minDist = Double.MAX_VALUE;
			int nextIndex = 0;
			for (Integer p:remainingChains){
				double distance = chainDist.get(lastIndex,p);
				if (distance<minDist){
					nextIndex = p;
					minDist = distance;
				}
			}
			//Add the chain atoms to the sorted array
			for (Atom a:chains.get(nextIndex)) sortedAtoms.add(a);
			remainingChains.remove((Integer) nextIndex);
			logger.info(chains.get(nextIndex)[0].getGroup().getChainId());
			lastIndex = nextIndex;
		}
		
		//Return the Atoms with the sorted chain order
		return sortedAtoms.toArray(new Atom[sortedAtoms.size()]);
	}
	
	/**
	 * Application: Cyclic Symmetry (Cn).<p>
	 * Assumes that the input Structure contains the chains 
	 * in the biological assembly, otherwise the asymmetric 
	 * unit symmetry will also be detected.
	 * 
	 * @param structure Structure containing the Chains
	 * @return sorted representative atom array
	 * @throws StructureException 
	 */
	public static Atom[] quatSort(Structure structure) 
			throws StructureException{
		
		QuatSymmetryParameters params = new QuatSymmetryParameters();
		QuatSymmetryDetector detector = 
				new QuatSymmetryDetector(structure, params);
		
		List<QuatSymmetryResults> global = detector.getGlobalSymmetry();
		QuatSymmetryResults result = global.get(0);
		
		if (result.getSymmetry().equals("C1")) { //asymmetric
			/*List<List<QuatSymmetryResults>> local = 
					detector.getLocalSymmetries();
			if (local.isEmpty())
				return StructureTools.getRepresentativeAtomArray(structure);
			else result = local.get(0).get(0);*/
			return StructureTools.getRepresentativeAtomArray(structure);
		} else if (!result.getSymmetry().contains("C")){ //non-cyclic
			return StructureTools.getRepresentativeAtomArray(structure);
		}
		
		List<Integer> chainOrder = new ArrayList<Integer>();
		int axisIndex = result.getRotationGroup().getPrincipalAxisIndex();
		List<Integer> perm = result.getRotationGroup().getRotation(axisIndex).
				getPermutation();
		int index = 0;
		
		//Follow the permutations to generate the cyclic order of the chains
		while (chainOrder.size() < perm.size()){
			chainOrder.add(index);
			index = perm.get(index);
		}
		
		List<Atom> atomList = new ArrayList<Atom>();
		for (int c : chainOrder){
			String id = result.getSubunits().getChainIds().get(c);
			Chain chain = structure.getChainByPDB(id);
			Atom[] array = StructureTools.getRepresentativeAtomArray(chain);
			for (Atom a : array) atomList.add(a);
		}
		
		return atomList.toArray(new Atom[atomList.size()]);
	}
	
	/**
	 * Test the chain ordering in some canonical examples:
	 * <ul><li>1KQ1: C6 Homo-hexamer, chain order A,B,M,K,I,H
	 * <li>4QVC: C6 Homo-hexamer, chain order A,B,C,D,E,F
	 * <li>4P24: C7 Homo-heptamer, chain order A,B,C,D,E,F,G
	 * <li>3X2R: C9 Homo, chain order A,B,C,D,E,F,G,H,I
	 * </ul>
	 */
	public static void main(String[] args) throws Exception {
		
		//String name = "4QVC";
		//String name  = "1KQ1";
		//String name = "4P24";
		//String name = "1uae";
		String name = "2VQA";
		
		//Load the biological assembly of the protein
		Structure structure = null;
		try {
			structure = StructureIO.getBiologicalAssembly(name);
		} catch (StructureException e) {
			structure = StructureIO.getBiologicalAssembly(name, 0);
		}
		
		Atom[] ca1 = ChainSorter.cyclicSort(structure);
		//Atom[] ca1 = ChainSorter.quatSort(structure);
		
		CeSymm cesymm = new CeSymm();
		CESymmParameters params = (CESymmParameters) cesymm.getParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setOptimization(true);
		params.setMultipleAxes(true);
		
		MultipleAlignment msa = cesymm.analyze(ca1);
		SymmetryDisplay.display(msa, cesymm.getSymmetryAxes());
	}
}
