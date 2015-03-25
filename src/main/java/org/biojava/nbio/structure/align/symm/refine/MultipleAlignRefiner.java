package org.biojava.nbio.structure.align.symm.refine;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Stack;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.subunit.SubunitTools;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;

/**
 * Creates a refined alignment with the multiple self-alignments obtained from blacking out previous alignments.
 * @author lafita
 */

public class MultipleAlignRefiner implements Refiner {

	public MultipleAlignRefiner() {
		super();
	}
	
	@Override
	public AFPChain refine(AFPChain[] afpAlignments, Atom[] ca1, Atom[] ca2, int order)
			throws RefinerFailedException {
		
		AFPChain refinedAFP = new AFPChain();
		
		try {
			 refinedAFP = cycleRefineAFP(afpAlignments, ca1);
		} catch (StructureException e) {
			e.printStackTrace();
		}
		
		return refinedAFP;
	}
	
	 /**
	 * Calculates from a set of AFP alignments the groups of residues (the group size is the order of symmetry)
	 * that align together and are consistent, equivalent for each subunit). As a result a modified AFPChain with 
	 * <order of symmetry> consistent groups is produced.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] array of the protein.
	 * OUTPUT: an AFP alignment with subunit groups consistent between each other.
	 * ALGORITHM: cycle detection for each residue checking each time that the distance between residues is sufficient
	 * RUNNING TIME: DFS of depth order, polynomial in length of the protein. No bottleneck.
	 */
	public static AFPChain cycleRefineAFP(AFPChain[] allAlignments, Atom[] ca1) throws StructureException {
		
		//Create the alignment graph and initialize a variable to store the groups
		//double[][][] weights = buildWeightedAFPgraph(allAlignments, ca1);
		List<List<Integer>> graph = buildAFPgraph(allAlignments, ca1);
		//Count the number of aligned residues to exclude the loops from the interresidue distances
		int aligned_res = 0;
		for (List<Integer> v:graph){
			//Count residues with order-1 edges (connected consistently in all the alignments)
			if (v.size()>0){
				aligned_res++;
			}
		}
		List<List<Integer>> groups = new ArrayList<List<Integer>>();
		
		//Variable to store the residues already present in one of the selected groups
		List<Integer> alreadySeen = new ArrayList<Integer>();
		int order = allAlignments.length+1;
		
		//Loop through all the residues in the graph (only consider the first residues, because cycles must include them to be consistent).
		for (int i=0; i<graph.size()/(order-order/2); i++){
			if (!alreadySeen.contains(i)){
			//System.out.println("Cycle for residue "+i);
			
			//Initialize the variables for the DFS of this residue iteration
			Stack<Integer> path = new Stack<Integer>(); //stack that stores the current path nodes
			Stack<List<Integer>> stack = new Stack<List<Integer>>(); //stack that stores the nodes to be visited next: [vertex,level]
			List<Integer> source = new ArrayList<Integer>(); //source information: level 0
			source.add(i);
			source.add(0);
			stack.push(source);
			
			boolean foundGroup = false; //Do not loop more if you already found a group of connected nodes
			
			while (!stack.isEmpty() && !foundGroup){
				
				/*//Print path and stack at each iteration
				System.out.println("Stack: ");
				for (List<Integer> s:stack){
					System.out.println(s);
				}
				
				System.out.println("Path: ");
				for (Integer s:path){
					System.out.println(s);
				}*/
				
				List<Integer> vertex = stack.pop();
				
				//If the vertex level is lower than the path size remove the last element of the path
				while (vertex.get(1)<=path.size()-1){
					//System.out.println("Popped from the path: "+path.pop());
					path.pop();
				}
				
				//If the vertex has level lower than the order consider its neighbors
				if (vertex.get(1)<order && !path.contains(vertex.get(0))){
					//First add the node to the path
					path.push(vertex.get(0));
					//System.out.println("Pushed to the path: "+path.peek());
					
					for (int k=0; k<graph.get(vertex.get(0)).size(); k++){
						//Extract the next node to be considered (neighbor k of the vertex)
						List<Integer> node = new ArrayList<Integer>();
						node.add(graph.get(vertex.get(0)).get(k));
						node.add(path.size());
						//Only add to the stack the nodes not included in the current path with level less than the order
						if (!path.contains(node.get(0)) && node.get(1)<order){
							stack.push(node);
						}
						//If the level=order and the node is equal to the source a cycle of size order has been found
						else if (node.get(0)==i && node.get(1)==order){
							//Initialize the group of residues
							List<Integer> group = new ArrayList<Integer>();
							int n = path.size();
							//Store the nodes in the path in the group and sort them
							for (int x=0; x<n; x++){
								int p = path.get(x);
								group.add(p);
							}
							Collections.sort(group);
							
							//Check that the residues have consistent interresidue distance (number of aligned residues in between) between the same group
							boolean consistent = true;
							
							for (int g=0; g<order; g++){
								for (int h=0; h<order; h++){
									int residueDist = Math.abs(group.get(g)-group.get(h));
									//Special case when comparing the first and last groups, because the length of the protein has to be considered
									if ((g==0 && h==order-1) || (h==0 && g==order-1)){
										residueDist = ca1.length - Math.abs(group.get(g)-group.get(h));
									}
									//If the distance between the two residues in number is lower they are too close in sequence and not consistent.
									if (residueDist<(aligned_res/(order+2)) && h!=g){
										consistent = false;
										//System.out.println("Not consistent: difference of "+residueDist+" with maximum "+aligned_res/(order+2)+"...");
										break;
									}
								}
							}
							
							//If any residue of the group is in another group already seen mark it as inconsistent
							for (int d=0; d<order; d++){
								if (alreadySeen.contains(group.get(d))){
									consistent = false;
								}
							}
							
							//Check that the group is consistent with the previous, all the residues should be greater than the last group
							int len = groups.size();
							if (len!=0){
								for (int d=0; d<order; d++){
									if (groups.get(len-1).get(d)>group.get(d)){
										//System.out.println("Inconsistent group: not increasing residues");
										consistent=false;
										break;
									}
								}
							}
							if (!consistent) continue;
							
							//If the conditions are fulfilled add the group
							groups.add(group);
							//System.out.println("Group added, size: "+group.size());
							for (int e:group){
								alreadySeen.add(e);
								//System.out.println(e);
							}
							foundGroup = true;
							path.clear();
							//System.out.println("Clearing path...");
							break;
						}
					}
				}
			} //end of DFS
			}
		} //end of all the residue analysis
		
		//Initialize the optAln variable
		List<List<List<Integer>>> optAln = new ArrayList<List<List<Integer>>>();
		for (int k=0; k<order; k++){
			List<List<Integer>> chains = new ArrayList<List<Integer>>();
			for (int j=0; j<2; j++){
				List<Integer> chain = new ArrayList<Integer>();
				chains.add(chain);
			}
			optAln.add(chains);
		}
		
		//Convert the groups of residues into the optimal alignment (suppose the groups are already sorted by their first residue)
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
		
		return AlignmentTools.replaceOptAln(optAlgn, allAlignments[order-2], ca1, ca1);
	}
	
	/**
	 * Calculates a graph in the format of adjacency list from the set of alignments, where each vertex is a 
	 * residue and each edge means the connection between the two residues in one of the alignments.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] array.
	 * OUTPUT: the alignment graph (describing relations between residues). List dimensions: AdjList[vertices][edges]
	 */
	public static List<List<Integer>> buildAFPgraph(AFPChain[] allAlignments, Atom[] ca1) {
		
		//Initialize the adjacency list that stores the graph
		List<List<Integer>> adjList = new ArrayList<List<Integer>>();
		for (int n=0; n<ca1.length; n++){
			List<Integer> edges = new ArrayList<Integer>();
			adjList.add(edges);
		}
		
		//Convert the multiple alignments into a triple list to handle the relations
		//List<List<List<Integer>>> alignments = processMultipleAFP(allAlignments);
		
		for (int k=0; k< allAlignments.length; k++){
			for (int i=0; i<allAlignments[k].getOptAln().length; i++){
				for (int j=0; j<allAlignments[k].getOptAln()[i][0].length; j++){
				
					//The vertex is the residue in the first chain and the edge the one in the second chain
					int vertex = allAlignments[k].getOptAln()[i][0][j];
					int edge = allAlignments[k].getOptAln()[i][1][j];
					if (!adjList.get(vertex).contains(edge)){
						adjList.get(vertex).add(edge);
					}
					//Make the graph undirected (optional feature)
					if (!adjList.get(edge).contains(vertex)){
						adjList.get(edge).add(vertex);
					}
				}
			}
		}
		/*//Print graph information
		System.out.println("GRAPH INFORMATION:");
		System.out.println(" - Vertices: "+adjList.size());
		int edges = 0;
		for (List<Integer> v:adjList){
			edges+=v.size();
		}
		System.out.println(" - Edges: "+edges);
		System.out.println(" - Consistent connected residues: "+count+" out of "+adjList.size());*/
		
		return adjList;
	}
	
	/**
	 * Displays the alignment of the protein with the refined optimal alignment with the subunit blocks.
	 * 
	 * Tried and worked for:
	 * 		- order 2: 2J5A, 3HDP, 4HHB, 1SQU.A, 2F9H.A, 3DDV.A, 4FI3.F, 1H9M.A, 1MP9.A, 1AIN, 1VYM.A, 4HPL.A, 1UBI, 1GUA.B, 1BA2.A, 1HIV
	 * 		- order 3: 4DOU, 1VYM, 2AFG.A, 1HCE, 1TIE, 4I4Q, 1itb.A, 3jut.A, d1jlya1
	 *   	- order 4: 1GEN, 1HXN
	 *  	- order 5: 1G61.A, 1TL2.A, 2jaj.A
	 *  	- order 6: 1U6D,
	 *      - order 7: 1JOF.A, 1JTD.B, 1K3I.A, 2I5I.A, 1JV2.A, 1GOT.B, 1A12.A
	 *      - order 8: 1TIM.A, 1VZW, 1NSJ
	 *      - helical: 1B3U.A, 1EZG.A, 1DFJ.I, 1AWC.B, 1D0B.A
	 *      - unknown: 1WD3, 1Z7X, 1DCE
	 *  	
	 * Did not work for:  1VYM.A (buggy rotation axis)
	 *                    3HKE.A (partial alignment)
	 *                    d1poqa_ (needs a threshold as a minimum score to proceed the multiple alignment, because the optimal is very low).
	 *                    
	 *                    
	 *                  
	 * BUGS:   1*- For the 1JTD.B structure, the blackout is not done properly in the last alignment and the alignment is made
	 *            over a black area, that is why the order of symmetry is incorrectly determined to 8 instead of 7. To reproduce
	 *            the error use CEsymmColorGUI, which displays the dotplot of the last alignment (same with 1TIM.A). Solved with the
	 *            matrix listener in the OrigM align method in the CeSymm class, that maintains the black regions during optimization.
	 *         2*- The 3D alignment deletes some regions in the rotated (second) protein, which may be the unaligned regions between
	 *            the subunits. It might be a problem with AlignmentTools.updateSuperposition. Some examples are: 1G61.A, 3HDP. Solved
	 *            by clonning the AFPChain instead of initializing it from 0, the method needs the original one.
	 *         3*- The information of gaps and RMSD in the Sequence Alignment Display is incorrect (it is set to 0). updateSuperposition 
	 *            has to be changed (the dummy code), in order to update the values.
	 *         4*- In the molecule 1G61.A the helices are not aligned although they seem to be symmetric, consider a less restrictive
	 *            last step to select them also as symmetric. Solved with a last step to select the groups that do not form cycles.
	 *         5*- The alignment panel does not color correctly the alignment, because it considers a " " in the symb alignment not as 
	 *            a mismatch (when the alignment is not FatCat only) but rather as an aligned pair and it colors them (modification in 
	 *            the biojava code DisplayAFP line 76), although they are not really aligned. Possible implications for biojava?
	 *         6*- When the subunits are colored in the 3D structure, in some structures the color does not change between two subunits,
	 *            they are either all blue or all green. This happens with the 1JTD.B structure (maybe red is not added properly).
	 *         7- Rotation axis is not displayed with an issue with the jmol window when evaluating the string. Examples: 3DDV.A
	 *         8*- The subunit selection seems to be very restrictive for proteins of higher order of symmetry. One solution could be
	 *            to consider, if there are not <order> cycles, cycles of smaller size and establish (or group) the subunits by pairwise
	 *            (or more) similarity groups. Different approach for order 6-8. Examples: 1VZW, 1TIM.A
	 *         9- For small proteins, the TM score is very low even for the first alignment, which results to incorrectly determine the 
	 *            order (higher than it is). This could be fixed by determining a threshold that considers the length of the protein.
	 *            Examples: 1GUA.B, 1UBI, d1poqa_.
	 *        10*- From the alignment panel, an error is thrown because a position cannot be matched. Example: 1A12.A. Unknown but solved.
	 *        11*- Protein 1G61.A gives some problems in the alignment: the colors are misplaced in the boundaries of the subunits and
	 *            the FATCAT result is not properly shown (it shows the text alignment instead). Something to do with a double gap that
	 *            might be ignore in the first subunit. Solved in getBlockNr of biojava display code.
	 *        12*- The getAlign method does not consider double gaps when calculating the alignment strings and that is why some errors
	 *            occur in the sequence alignment Display (color not correct and repeated residues). The problem was actually in the
	 *            OptAln, because the residues were not contiguous in all the subunits. Solved by checking consistency between groups in
	 *            the refinement method, there was a bug in the names of variables.
	 *        13- In 1VYM structure there is a loop identified as not aligned, but it is present in the three subunits and the sequence
	 *            is highly conserved, although in the 3D alignment is only seems to align well two of the three loops. In this case the
	 *            subunit conditions are too restrictive.
	 *         
	 *         * Solved!
	 *                    
	 *  
	 * @author lafita
	 *
	 */
	public static void main(String[] args) throws StructureException, IOException{
		
		//String[] names = {"2F9H.A", "1SQU.A", "3HDP", "2AFG.A", "4DOU", "1GEN", "1G61.A", "1U6D", "1JOF.A", "1JTD.B", "1TL2.A", "2I5I.A", "1GOT.B", "1VZW", "1NSJ");
		//String[] names = {"d1poqa_", "1itb.A", "3jut.A", "2jaj.A", "d1jlya1" ,"1hiv"}; //New structures to test
		String[] names = {"d1poqa_"};
		
		for (int i=0; i<names.length; i++){
			
			//Set the name of the protein structure to analyze
			System.out.println("Analyzing protein "+names[i]);
			AtomCache cache = new AtomCache();
			String name = names[i];

			//Parse atoms of the protein into two DS
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			System.out.println("Protein length: "+ca1.length);
			
			//Initialize a new CeSymm class and its parameters and a new alignment class
			CeSymm ceSymm = new CeSymm();
			AFPChain afpChain = new AFPChain();
			
			//Perform the alignment and store
			afpChain = ceSymm.align(ca1, ca2);
			
			afpChain.setName1(name);
			afpChain.setName2(name);
				
			//Display the AFP alignment of the subunits
			StructureAlignmentJmol jmolPanel;
			jmolPanel = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
				
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(afpChain);
			jmolPanel.evalString(axis.getJmolScript(ca1));
			
			//Display two subunits superimposed
			//SubunitTools.displaySuperimposedSubunits(afpChain, name, ca1, ca2);
		}
	}

}
