package org.biojava.nbio.structure.align.symm;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * This class provides methods for sorting the chains of a structure in the symmetry order,
 * so that the CeSymm algorithm can deal with quaternary symmetry.
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
			chainAtoms.add(atoms);
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
	 * In case of homo n-meric cyclic symmetry (Cn) the following theorem applies: 
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
			if (debug) System.out.println(chains.get(nextIndex)[0].getGroup().getChainId());
			lastIndex = nextIndex;
		}
		
		//Return the Atoms with the sorted chain order
		return sortedAtoms.toArray(new Atom[sortedAtoms.size()]);
	}
	
	/**
	 * Test the chain ordering in some canonical examples:
	 * <ul><li>4PJO: Hetero-heptamer, chain order B,A,G,E,F,D,C
	 * <li>3CW1: Hetero-heptamer, chain order 
	 * <li>1KQ1: Homo-hexamer, chain order A,B,M,K,I,H
	 * </ul>
	 */
	public static void main(String[] args) throws Exception {
		
		//String name  = "3CW1.A:,B:,C:,D:,E:,F:,G:";
		//String name = "4PJO.A:,B:,C:,D:,E:,F:,G:";
		String name  = "1KQ1.A:,B:,H:,I:,K:,M:";
		
		AtomCache cache = new AtomCache();
		Structure structure = cache.getStructure(name);
		
		Atom[] ca1 = cyclicSorter(structure);
		//Atom[] ca1 = StructureTools.getRepresentativeAtomArray(structure);
		Atom[] ca2 = StructureTools.cloneAtomArray(ca1);
		
		CeSymm cesymm = new CeSymm();
		CESymmParameters params = (CESymmParameters) cesymm.getParameters();
		//params.setRefineMethod(RefineMethod.SINGLE);
		//params.setOptimization(true);
		
		AFPChain afp = cesymm.align(ca1, ca2);
		
		StructureAlignmentDisplay.display(afp, ca1, ca2);
		
	}
}
