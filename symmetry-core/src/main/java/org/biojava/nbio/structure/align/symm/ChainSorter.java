package org.biojava.nbio.structure.align.symm;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.symmetry.gui.SymmetryDisplay;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;

/**
 * This class provides methods for sorting the chains 
 * of a structure in their symmetric order, so that 
 * the CeSymm algorithm (and other sequence-dependent 
 * algorithms) can deal with multiple chains.
 * 
 * @author Aleix Lafita
 *
 */
public class ChainSorter {

	private static final boolean debug = true;
	
	/**
	 * Application: Cyclic Symmetry (Cn).
	 * <p>
	 * This method sorts the chains by the following steps:
	 * <ul><li>Calculate the centroid of each chain.
	 * <li>Pick the farthest chain from all others as the first (for open cases, can be random for closed).
	 * <li>Iteratively pick the closest chain to the last one from the remaining set.
	 * </ul>
	 * In case of homo n-meric cyclic symmetry (Cn) the following theorem applies: 
	 * "the closest point of another point is either the previous or the next in the 
	 * 3D rotation around the axis of symmetry". That is true because the points are situated at
	 * the edge of a circle. Thus, this property can be used to sort the chains in cyclic order.
	 * 
	 * @param structure Structure containing the Chains
	 * @return Atom[] with the sorted order of Atoms, corresponding to the ordering of the chains.
	 */
	public static Atom[] cyclicSorter(Structure structure){
		
		List<Atom[]> chainAtoms = new ArrayList<Atom[]>();
		for (Chain c:structure.getChains()){
			Atom[] atoms = StructureTools.getRepresentativeAtomArray(c);
			if (atoms.length > 0) chainAtoms.add(atoms);
		}
		return cyclicSorter(chainAtoms);
	}
	
	/**
	 * Application: Cyclic Symmetry (Cn).
	 * <p>
	 * This method sorts the chains by the following steps:
	 * <ul><li>Calculate the centroid of each chain.
	 * <li>Pick the farthest chain from all others as the first (for open cases, can be random for closed).
	 * <li>Iteratively pick the closest chain to the last one from the remaining set.
	 * </ul>
	 * In case of n-meric cyclic symmetry (Cn) the following theorem applies: 
	 * "the closest point of another point is either the previous or the next in the 
	 * 3D rotation around the axis of symmetry". That is true because the points are situated at
	 * the edge of a circle. Thus, this property can be used to sort the chains in cyclic order.
	 * 
	 * @param chains List of Atoms of each chain
	 * @return Atom[] with the sorted order of Atoms, corresponding to the ordering of the chains.
	 */
	public static Atom[] cyclicSorter(List<Atom[]> chains){
		
		if (debug) System.out.println("Ordering of "+chains.size()+" chains:");
		//Initialize the variables nedeed
		List<Point3d> centroids = new ArrayList<Point3d>();
		Matrix chainDist = new Matrix(chains.size(), chains.size());
		List<Atom> sortedAtoms = new ArrayList<Atom>();
		List<Integer> remainingChains = new ArrayList<Integer>();
		for (int n=0; n<chains.size(); n++) remainingChains.add(n); //Contains the indices of the chains not seen
		
		/** STEP 1: Centroid calculation */
		//Loop thorugh all the Atom arrays of the chains and calculate the centroids
		for (Atom[] array:chains){
			Point3d centroid = new Point3d();
			//Loop through all atoms in the chain and get the average distance
			for (Atom a:array){
				centroid.x += a.getX();
				centroid.y += a.getY();
				centroid.z += a.getZ();
			}
			centroid.x /= array.length;
			centroid.y /= array.length;
			centroid.z /= array.length;
			centroids.add(centroid);
		}
		
		/** STEP 2: Pick the farthest chain to all others */
		//Calculate the interchain chain distances
		for (int m=0; m<centroids.size(); m++){
			for (int n=0; n<centroids.size(); n++){
				chainDist.set(m, n, centroids.get(m).distance(centroids.get(n)));
			}
		}
		//Pick the farthest point to all others
		double maxDist = 0.0;
		int maxIndex = 0;
		for (int p=0; p<centroids.size(); p++){
			double distance = 0.0;
			for (int p2=0; p2<centroids.size(); p2++) distance += chainDist.get(p, p2);
			if (distance>maxDist){
				maxIndex = p;
				maxDist = distance;
			}
		}
		
		/** STEP 3: Iteratively pick the closest chain to the last one */
		//First add the farthest chain to the beginning of the array
		for (Atom a:chains.get(maxIndex)) sortedAtoms.add(a);
		remainingChains.remove((Integer) maxIndex);
		if (debug) System.out.println(chains.get(maxIndex)[0].getGroup().getChainId());
		
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
			if (debug) 
				System.out.println(chains.get(nextIndex)[0].getGroup().getChainId());
			lastIndex = nextIndex;
		}
		
		//Return the Atoms with the sorted chain order
		return sortedAtoms.toArray(new Atom[sortedAtoms.size()]);
	}
	
	/**
	 * Application: any quaternary symmetry supported.<p>
	 * Assumes that the input Structure contains the chains 
	 * in the biological assembly, otherwise the asymmetric 
	 * unit symmetry will also be detected.
	 * 
	 * @param structure Structure containing the Chains
	 * @return
	 */
	public static Atom[] quaternaryAxisSorter(Structure structure){
		
		return null;
	}
	
	/**
	 * Test the chain ordering in some canonical examples:
	 * <ul><li>4PJO: C7 Hetero-heptamer, chain order B,A,G,E,F,D,C
	 * <li>3CW1: C7 Hetero-heptamer, chain order C,F,E,G,D,A,B
	 * <li>1KQ1: C6 Homo-hexamer, chain order A,B,M,K,I,H
	 * <li>2WBH: A180 Icosahedral
	 * <li>1AEW: A24 octahedral
	 * </ul>
	 */
	public static void main(String[] args) throws Exception {
		
		//String name  = "3CW1";
		//String name = "4PJO";
		String name  = "1KQ1";
		//String name = "1AEW";
		
		//Load the biological assembly of the protein
		Structure structure = null;
		try {
			structure = StructureIO.getBiologicalAssembly(name, 1);
		} catch (StructureException e) {
			structure = StructureIO.getBiologicalAssembly(name, 0);
		}
		
		Atom[] ca1 = cyclicSorter(structure);
		
		CeSymm cesymm = new CeSymm();
		CESymmParameters params = (CESymmParameters) cesymm.getParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setOptimization(true);
		
		MultipleAlignment msa = cesymm.analyze(ca1);
		SymmetryDisplay.display(msa);
	}
}
