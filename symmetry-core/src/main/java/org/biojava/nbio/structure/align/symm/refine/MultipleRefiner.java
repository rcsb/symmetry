package org.biojava.nbio.structure.align.symm.refine;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.symmetry.internal.Refiner;
import org.biojava.nbio.structure.symmetry.internal.RefinerFailedException;

/**
 * Creates a refined alignment with the multiple self-alignments 
 * obtained from blacking out previous alignments.
 * Still in development, because the clustering is not robust enough.
 * 
 * @author Aleix Lafita
 * 
 */
public class MultipleRefiner implements Refiner {

	@Override
	public AFPChain refine(List<AFPChain> afpAlignments, Atom[] atoms,
			int order) throws StructureException, RefinerFailedException {

		//return cycleRefine(afpAlignments, atoms, order);
		throw new RefinerFailedException("Multiple Refiner not yet implemented!");
	}

	/**
	 * Calculates from a set of AFP alignments the groups of residues 
	 * (the group size is the order of symmetry) that align together 
	 * and are consistent, equivalent for each subunit).
	 * <p>
	 * As a result a modified AFPChain with  <order of symmetry> 
	 * consistent groups is returned.
	 * <p>
	 * It uses a DFS into the symmetry graph to find clusters of
	 * consistently aligned residues.
	 */
	private static AFPChain cycleRefine(List<AFPChain> allAlignments, 
			Atom[] atoms, int order) throws StructureException {

		//TODO implement again when jgrapht added as dependency
		/*Graph<Integer> graph = 
				SymmetryTools.buildSymmetryGraph(allAlignments, atoms);
		List<List<Integer>> groups = new ArrayList<List<Integer>>();

		//Variable to store the residues already present in one of the groups
		List<Integer> alreadySeen = new ArrayList<Integer>();

		for (int i=0; i<graph.size(); i++){
			if (!alreadySeen.contains(i)){

			//Initialize the variables for the DFS of this residue iteration
			Stack<Integer> path = new Stack<Integer>();
			Stack<List<Integer>> stack = new Stack<List<Integer>>();
			List<Integer> source = new ArrayList<Integer>();
			source.add(i);
			source.add(0);
			stack.push(source);

			boolean foundGroup = false;

			while (!stack.isEmpty() && !foundGroup){

				List<Integer> vertex = stack.pop();

				while (vertex.get(1)<=path.size()-1){
					path.pop();
				}

				//consider its neighbors
				if (vertex.get(1)<order && !path.contains(vertex.get(0))){
					//First add the node to the path
					path.push(vertex.get(0));
					List<Integer> neighbors = 
							graph.getNeighborIndices(vertex.get(0));

					for (int k=neighbors.size()-1; k>=0; k--){
						//Extract the next node to be considered
						List<Integer> node = new ArrayList<Integer>();
						node.add(neighbors.get(k));
						node.add(path.size());
						if (!path.contains(node.get(0)) 
								&& node.get(1) < order 
								&& !alreadySeen.contains(node.get(0))){
							stack.push(node);
						}//cycle of size order has been found
						else if (node.get(0)==i && node.get(1)==order){
							//Initialize the group of residues
							List<Integer> group = new ArrayList<Integer>();
							int n = path.size();
							//Store the nodes in the path in the group
							for (int x=0; x<n; x++){
								int p = path.get(x);
								group.add(p);
							}
							Collections.sort(group);

							boolean consistent = true;
							for (int g=0; g<order; g++){
								for (int h=0; h<order; h++){
									int d = Math.abs(group.get(g)
											-group.get(h));
									if ((g==0 && h==order-1) 
											|| (h==0 && g==order-1)){
										d = atoms.length - 
												Math.abs(group.get(g)-
														group.get(h));
									}
									if (d<(atoms.length/(order+2)) && h!=g){
										consistent = false;
										break;
									}
								}
							}

							int len = groups.size();
							if (len!=0){
								for (int d=0; d<order; d++){
									if (groups.get(len-1).get(d)>group.get(d)){
										consistent=false;
										break;
									}
								}
							}
							if (!consistent) continue;

							//If the conditions are fulfilled add the group
							groups.add(group);
							for (int e:group){
								alreadySeen.add(e);
							}
							foundGroup = true;
							path.clear();
							break;
						}
					}
				}
			} //end of DFS
			}
		} //end of all the residue analysis

		//Initialize the optAln variable
		List<List<List<Integer>>> optAln = 
				new ArrayList<List<List<Integer>>>();
		for (int k=0; k<order; k++){
			List<List<Integer>> chains = new ArrayList<List<Integer>>();
			for (int j=0; j<2; j++){
				List<Integer> chain = new ArrayList<Integer>();
				chains.add(chain);
			}
			optAln.add(chains);
		}

		//Convert the groups of residues into the optimal alignment
		for (List<Integer> group:groups){
			//System.out.println("Group: ");
			for (int k=0; k<group.size(); k++){
				optAln.get(k).get(0).add(group.get(k));
				optAln.get(k).get(1).add(group.get((k+1)%order));
				//System.out.println(group.get(k));
			}
		}

		//Convert the triple list into a triple array of the AFPChain format
		int[][][] optAlgn = new int[order][2][];
		for (int i=0; i<order; i++){
			for (int k=0; k<2; k++){
				int n = optAln.get(i).get(0).size();
				optAlgn[i][k] = new int[n];
				for (int j=0; j<n; j++){
					optAlgn[i][k][j] = optAln.get(i).get(k).get(j);
				}
			}
		}

		return AlignmentTools.replaceOptAln(
				optAlgn, allAlignments.get(order-2), atoms, atoms);*/
		return null;
	}
}
