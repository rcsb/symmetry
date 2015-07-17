package org.biojava.nbio.structure.align.symm.refine;

import java.util.List;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.order.OrderDetector;

/**
 * Creates a refined alignment with the multiple self-alignments 
 * obtained from blacking out previous alignments.
 * Still in development, because the clustering is not robust enough.
 * 
 * @author Aleix Lafita
 * 
 */
public class MultipleRefiner implements Refiner {

	private OrderDetector orderDetector;

	/**
	 * The class needs an orderDetector because the order might differ 
	 * among the different alignments.
	 * 
	 * @param orderDetector
	 */
	public MultipleRefiner(OrderDetector orderDetector) {
		this.orderDetector = orderDetector;
	}

	@Override
	public AFPChain refine(List<AFPChain> afpAlignments, Atom[] atoms, 
			int order) throws RefinerFailedException,StructureException {

		//Use the help of the SINGLE refiner to increase the consistency of the alignments
		for (int i=0; i<afpAlignments.size(); i++){
			try {
				order = orderDetector.calculateOrder(afpAlignments.get(i), atoms);
				afpAlignments.set(i,SingleRefiner.refineSymmetry(afpAlignments.get(i), atoms, atoms, order));
			} //It means that one of the refined alignments has length 0, so continue without refining
			catch (Exception ignore){
				//afpAlignments[i] = null;
				continue;
			}
		}
		//return cycleRefine(afpAlignments, ca1, ca2, order);
		return maxOrder(afpAlignments);
	}

	/**
	 *  Returns the alignment with the maximum length of the set of afpAlignments.
	 */
	private AFPChain maxLength(List<AFPChain> afpAlignments) {

		AFPChain maxAFP = null;
		for (int i=0; i<afpAlignments.size(); i++){
			if (afpAlignments.get(i)!=null){
				if (maxAFP==null) maxAFP = afpAlignments.get(i);
				else if (maxAFP.getOptLength()<afpAlignments.get(i).getOptLength()) maxAFP = afpAlignments.get(i);
			}
		}
		return maxAFP;
	}

	/**
	 *  Returns the alignment with the maximum TM score of the set of afpAlignments.
	 */
	private AFPChain maxScore(List<AFPChain> afpAlignments) {

		AFPChain maxAFP = null;
		for (int i=0; i<afpAlignments.size(); i++){
			if (afpAlignments.get(i)!=null){
				if (maxAFP==null) maxAFP = afpAlignments.get(i);
				else if (maxAFP.getTMScore()<afpAlignments.get(i).getTMScore()) maxAFP = afpAlignments.get(i);
			}
		}
		return maxAFP;
	}

	/**
	 *  Returns the signifcant alignment with the maximum order of the set of afpAlignments (and maximum TM score
	 *  if there are more than one with the same order).
	 */
	private AFPChain maxOrder(List<AFPChain> afpAlignments) {

		AFPChain maxAFP = null;
		for (int i=0; i<afpAlignments.size(); i++){
			if (afpAlignments.get(i)!=null){
				if (maxAFP==null) maxAFP = afpAlignments.get(i);
				else if (maxAFP.getBlockNum()<afpAlignments.get(i).getBlockNum() && afpAlignments.get(i).getTMScore() > CeSymm.symmetryThreshold-0.1)
					maxAFP = afpAlignments.get(i);
				else if (maxAFP.getBlockNum()==afpAlignments.get(i).getBlockNum())
					if (maxAFP.getTMScore()<afpAlignments.get(i).getTMScore()) maxAFP = afpAlignments.get(i);
			}
		}
		return maxAFP;
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
