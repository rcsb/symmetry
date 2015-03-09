package org.biojava.nbio.structure.align.symm.subunit;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.align.util.AlignmentTools;

/**
 * Methods to process AFPChain multiple alignments and analyze the symmetrical subunits generated.
 * 
 * @author Aleix Lafita
 * 
 * Last modified: 09.03.2015
 */
public class SubunitTools {
	
	/**
	 * NOT USED ANYMORE: new method refinedAFP outputs an AFP with the subunit residues as blocks.
	 * 
	 * Method that converts a set of alignments (corresponding to the symmetry rotations) into a list of residues 
	 * for each subunit of the protein that is consistent (equivalent starting and ending residues for subunit).
	 * 
	 * INPUT: a list of AFPChain alignments and the protein atoms.
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It extracts the subunit intervals from the alignments and then includes as part of every subunit only the 
	 *            residues that are present in all the alignment and are consistent.
	 */
	public static List<List<Integer>> extractSubunits(Atom[] ca1, List<AFPChain> allAlignments){
	
		//Process the alignments to which the subunit residues should be restricted (number higher is more restrictive)
		List<List<Integer>> afpResidues = new ArrayList<List<Integer>>();
		int number = allAlignments.size(); //can be in the interval [1,allAlignments.size()]
		for (int k=0; k<number; k++){
			List<Integer> residues = residuesAFP(allAlignments.get(k));
			afpResidues.add(residues);
		}
				
		//Calculate the intervals
		List<Integer> intervals = calculateIntervals(ca1, allAlignments);
		
		//Define the residues of each subunit with the intervals (cut regions for the subunits) and only add each residue
		//if it is contained in all the alignments
		List<List<Integer>> subunits = new ArrayList<List<Integer>>();
		List<Integer> subunit = new ArrayList<Integer>();
		//sub variable keeps track of the current subunit
		int sub = 1;
		for (int i=0; i<ca1.length; i++){
						
			//If the residue is present in all the alignments include it to the subunit
			boolean include = true;
			for (int k=0; k<afpResidues.size();k++){
				if (!afpResidues.get(k).contains(i)){
					include = false;
					break;
				}
			}
			if (include){
				subunit.add(i);
			}
			
			//If the residue is one of the subunit cutting points store its residues and go to the next subunit
			if (intervals.get(sub)==i){
				List<Integer> copySub = new ArrayList<Integer>(subunit);
				subunits.add(copySub);
				subunit.clear();
				if (sub==intervals.size()-1)
					break;
				i=intervals.get(sub+1);
				sub+=2;
			}
		}
		return subunits;
	}

	/**
	 * Method that calculates the intervals of every subunit (cutting regions). 
	 * 
	 * INPUT: a list of AFP alignments and the array of atoms of the protein.
	 * OUTPUT: a list containing the cutting intervals (size: order+1).
	 * ALGORITHM: It takes all the starting and ending residues from the alignment blocks and clusters them by iteratively
	 *            deleting the higher of the closest pair of residues. Other approaches to cluster are being considered.
	 *            Need to be improved in the future, not using the HEURISTIC of most distant pairs.
	 * 
	 */
	public static List<Integer> calculateIntervals(Atom[] ca1, List<AFPChain> allAlignments){
				
		//Overlapping intervals problem: we have a set of 4*n intervals that overlap and define n subunits
		//Try to define the n subunits from the residue sets by its starting and ending residues of each interval
		List<Integer> intervals = new ArrayList<Integer>();
		int order = allAlignments.size()+1;
		
		//Add the first and last residue numbers of the protein
		intervals.add(0);
		intervals.add(ca1.length-1);
		
		for (int k=0; k<(order-1); k++){
			for (int y=0; y<allAlignments.get(k).getOptAln().length; y++){
				for (int z=0; z<allAlignments.get(k).getOptAln()[y].length; z++){
					//Add the first and last residue numbers of each alignment block
					int n = allAlignments.get(k).getOptAln()[y][z].length;
					if (!intervals.contains(allAlignments.get(k).getOptAln()[y][z][0]))
						intervals.add(allAlignments.get(k).getOptAln()[y][z][0]);
					if (!intervals.contains(allAlignments.get(k).getOptAln()[y][z][n-1]))
						intervals.add(allAlignments.get(k).getOptAln()[y][z][n-1]);
				}
			}
		}
		
		//Sort the interval numbers and take the most distant ones until there are <2*order> numbers
		intervals.sort(null);
		List<Integer> selectedInter = new ArrayList<Integer>();
		List<Integer> alreadySeen = new ArrayList<Integer>();
		while (selectedInter.size()<2*order){
			int selected = 1;
			int maxDistance = -1;
			int distance = 0;
			for (int j=1; j<intervals.size(); j++){
				distance = intervals.get(j)-intervals.get(j-1);
				if (distance > maxDistance && !alreadySeen.contains(j)){
					maxDistance = distance;
					selected = j;
				}
			}
			alreadySeen.add(selected);
			selectedInter.add(intervals.get(selected));
			selectedInter.add(intervals.get(selected-1));
		}
		selectedInter.sort(null);
		
		//Temporal: analyze the intervals
		System.out.println("Number of interval residues: "+selectedInter.size());
		for (int intr:selectedInter){
			System.out.println(intr);
		}
		
		return selectedInter;
	}
	
	/**
	 * NOT USED ANYMORE: replaced like the extractSubunits method.
	 * 
	 * Method that returns a list of the residues present in both chains of the AFP alignment.
	 * 
	 * INPUT: an AFP alignment.
	 * OUTPUT: a list containing all the residues in the AFP alignment, sorted increasingly.
	 */
	public static List<Integer> residuesAFP(AFPChain afpChain){
		
		//Get all the residue numbers in the alignment and store them in a list
		List<Integer> residues = new ArrayList<Integer>();
		for (int y=0; y<afpChain.getOptAln().length; y++){
			
			//Get the total number of residues in the group
			int n = afpChain.getOptAln()[y][0].length;
			
			//Append every residue to a list of residues
			for (int j=0; j<n; j++){
				residues.add(afpChain.getOptAln()[y][0][j]);
			}
		}
		
		//Now do the same for the other chain, select the residues only if present in the two chains.
		List<Integer> intersectionResidues = new ArrayList<Integer>();
		for (int y=0; y<afpChain.getOptAln().length; y++){
			
			//Get the total number of residues in the group
			int n = afpChain.getOptAln()[y][1].length;
			
			//Append every residue to a list of residues only if contained in the other rotation as well
			for (int j=0; j<n; j++){
				int residue = afpChain.getOptAln()[y][1][j];
				if (residues.contains(residue)){
					intersectionResidues.add(residue);
				}
			}
		}
		//Sort the residues increasingly
		intersectionResidues.sort(null);
		return intersectionResidues;
	}
	
	/**
	 * Method that processes an AFP alignment and returns a triple list of aligned residues.
	 * 
	 * INPUT: an list AFP alignment.
	 * OUTPUT: a triple list: first index = alignment number, second = chain number, third = position in the alignment.
	 */
	public static List<List<List<Integer>>> processMultipleAFP(List<int[][][]> allAlignments){
		
		//Initialize the triple list to be returned
		List<List<List<Integer>>> alignments = new ArrayList<List<List<Integer>>>();
		
		//Loop through all the alignments and extract a double list of residue pairs
		for (int k=0; k<allAlignments.size(); k++){
			
			//This list contains a list for each chain of the alignments with its residues
			List<List<Integer>> alignment = new ArrayList<List<Integer>>();
			
			for (int x=0; x<2; x++){
				//Get all the residues in one chain of the alignment and add them to the alignment list
				ArrayList<Integer> residues = new ArrayList<Integer>();
				for (int y=0; y<allAlignments.get(k).length; y++){
					
					//Get the total number of residues in the group
					int n = allAlignments.get(k)[y][x].length;
					//Append every residue to a list of residues
					for (int j=0; j<n; j++){
						residues.add(allAlignments.get(k)[y][x][j]);
					}
				}
				alignment.add(residues);
			}
			alignments.add(alignment);
		}
		return alignments;
	}
	
	/**
	 * Method that displays two or all of the superimposed subunits in jmol. In construction...
	 * 
	 * INPUT: a list of subunits and a list of AFP alignments.
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It extracts the subunit intervals from the alignments and then includes as part of every subunit only the 
	 *            residues that are present in all the alignments that contain that subunit.
	 */
	public static void displaySuperimposedSubunits(ArrayList<ArrayList<Integer>> subunits, ArrayList<AFPChain> allAlignments){
		
		//TODO
	}
	
	/**
	 * Method that displays the optimal structural alignment and colors each subunit differently in jmol. In construction...
	 * 
	 * INPUT: a list of AFP
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It displays the protein with jmol and colors the subunit residues using jmol commands.
	 */
	public static void displayColorSubunits(ArrayList<ArrayList<Integer>> subunits, String name){
		
		//TODO
	}
	
	/**
	 * It returns an AFPChain with a segmented optimal alignment, which means that it has <order of symmetry> blocks of
	 * aligned subunits.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] arrays.
	 * OUTPUT: the optimal AFPChain alignment divided into the subunits.
	 */
	public static AFPChain createOptAln(ArrayList<AFPChain> allAlignments, Atom[] ca1) throws StructureException {
		
		//Extract the subunits from the set of alignments
		List<List<Integer>> subunits = extractSubunits(ca1, allAlignments);
		int order = subunits.size();
		
		AFPChain optimalAFP = allAlignments.get(0);
				
		//Extract the old alignment
		int[][][] oldAlgn = optimalAFP.getOptAln();
		
		//Create a triple list with the optimal alignment divided into the subunits in <order> blocs
		ArrayList<ArrayList<ArrayList<Integer>>> newAlgn = new ArrayList<ArrayList<ArrayList<Integer>>>();
		
		//Loop through all the optimal alignment and subunits and select the pairs of aligned residues for each subunit
		for (int i=0; i<order; i++){
			ArrayList<Integer> chain1 = new ArrayList<Integer>();
			ArrayList<Integer> chain2 = new ArrayList<Integer>();
			for (int k=0; k<oldAlgn.length; k++){
				for (int j=0; j<oldAlgn[k][0].length; j++){
					if (subunits.get(i).contains(oldAlgn[k][0][j])){
						chain1.add(oldAlgn[k][0][j]);
						chain2.add(oldAlgn[k][1][j]);
					}
				}
			}
			ArrayList<ArrayList<Integer>> subunit = new ArrayList<ArrayList<Integer>>();
			subunit.add(chain1);
			subunit.add(chain2);
			newAlgn.add(subunit);
		}
				
		//Calculate the alignment length from all the subunits lengths
		int[] optLens = new int[order];
		for(int s=0;s<order;s++) {
			optLens[s] = newAlgn.get(s).get(0).size();
		}
		int optLength = 0;
		for(int s=0;s<order;s++) {
			optLength += optLens[s];
		}
		
		//Convert the three dimensional list newAlgn into a three dimensional array
		int[][][] optAlgn = new int[order][2][];
		for (int i=0; i<order; i++){
			for (int k=0; k<2; k++){
				optAlgn[i][k] = new int[optLens[i]];
				for (int j=0; j<optLens[i]; j++){
					optAlgn[i][k][j] = newAlgn.get(i).get(k).get(j);
				}
			}
		}
		
		//Temporal: print the sizes to check correctness
		System.out.println("Number of subunits: "+optAlgn.length);
		for (int i=0; i<optAlgn.length; i++){
			System.out.println("Subunit length: "+optAlgn[i][0].length);
		}
		
		//Create a new AFPChain and set everything needed from the optimal alignment (allAlignments.get(0))
		AFPChain refinedAFP = new AFPChain();
		refinedAFP.setCa1Length(optimalAFP.getCa1Length());
		refinedAFP.setCa2Length(optimalAFP.getCa2Length());
		
		//Set the new parameters of the optimal alignment
		refinedAFP.setOptLength(optLength);
		refinedAFP.setOptLen(optLens);
		refinedAFP.setOptAln(optAlgn);
		refinedAFP.setBlockNum(order);
		
		//Recalculate properties: superposition, tm-score, etc
		//Atom[] ca1clone = StructureTools.cloneCAArray(ca1); // don't modify ca1 positions
		AlignmentTools.updateSuperposition(refinedAFP, ca1, ca1);
		
		//It re-does the sequence alignment strings from the OptAlgn information only (why does not work..?)
		AFPAlignmentDisplay.getAlign(refinedAFP, ca1, ca1);
				
		return refinedAFP;
	}
	
	/**
	 * It returns an AFPChain with a segmented optimal alignment, which means that it has <order of symmetry> blocks of
	 * aligned subunits.
	 * 
	 * INPUT: the optimal alignment in a triple list (same format as optAln of AFPChain) and the Atom[] arrays.
	 * OUTPUT: the optimal AFPChain alignment object divided into the subunits.
	 */
	public static AFPChain createOptAln(List<List<List<Integer>>> newAlgn, Atom[] ca1) throws StructureException {
		
		//The order is the number of groups in the newAlgn
		int order = newAlgn.size();
		
		//Calculate the alignment length from all the subunits lengths
		int[] optLens = new int[order];
		for(int s=0;s<order;s++) {
			optLens[s] = newAlgn.get(s).get(0).size();
		}
		int optLength = 0;
		for(int s=0;s<order;s++) {
			optLength += optLens[s];
		}
		
		//Convert the three dimensional list newAlgn into a three dimensional array
		int[][][] optAlgn = new int[order][2][];
		for (int i=0; i<order; i++){
			for (int k=0; k<2; k++){
				optAlgn[i][k] = new int[optLens[i]];
				for (int j=0; j<optLens[i]; j++){
					optAlgn[i][k][j] = newAlgn.get(i).get(k).get(j);
				}
			}
		}
		
		//Temporal: print the sizes to check correctness
		System.out.println("Number of subunits: "+optAlgn.length);
		for (int i=0; i<optAlgn.length; i++){
			System.out.println("Subunit length: "+optAlgn[i][0].length);
		}
		
		//Create a new AFPChain and set everything needed from the optimal alignment (allAlignments.get(0))
		AFPChain refinedAFP = new AFPChain();
		refinedAFP.setCa1Length(ca1.length);
		refinedAFP.setCa2Length(ca1.length);
		
		//Set the new parameters of the optimal alignment
		refinedAFP.setOptLength(optLength);
		refinedAFP.setOptLen(optLens);
		refinedAFP.setOptAln(optAlgn);
		refinedAFP.setBlockNum(order);
		
		//Recalculate properties: superposition, tm-score, etc
		//Atom[] ca1clone = StructureTools.cloneCAArray(ca1); // don't modify ca1 positions
		AlignmentTools.updateSuperposition(refinedAFP, ca1, ca1);
		
		//It re-does the sequence alignment strings from the OptAlgn information only (why does not work..?)
		AFPAlignmentDisplay.getAlign(refinedAFP, ca1, ca1);
				
		return refinedAFP;
	}
	
	/**
	 * Calculates a directed graph in the format of adjacency list from the set of alignments, where each vertex is a 
	 * residue and each edge means the connection between the two residues in one of the alignments.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] array.
	 * OUTPUT: the alignment graph (describing relations between residues). List dimensions: AdjList[vertices][edges]
	 */
	public static List<List<Integer>> buildAFPgraph(List<int[][][]> allAlignments, Atom[] ca1) {
		
		//Initialize the adjacency list that stores the graph
		List<List<Integer>> adjList = new ArrayList<List<Integer>>();
		for (int n=0; n<ca1.length; n++){
			List<Integer> edges = new ArrayList<Integer>();
			adjList.add(edges);
		}
		
		//Convert the multiple alignments into a triple list to handle the relations
		List<List<List<Integer>>> alignments = processMultipleAFP(allAlignments);
		
		for (int k=0; k< alignments.size(); k++){
			for (int i=0; i<alignments.get(k).get(0).size(); i++){
				
				//The vertex is the residue in the first chain and the edge the one in the second chain
				int vertex = alignments.get(k).get(0).get(i);
				int edge = alignments.get(k).get(1).get(i);
				if (!adjList.get(vertex).contains(edge)){
					adjList.get(vertex).add(edge);
				}
			}
		}
		return adjList;
	}
	
	/**
	 * Calculates from a set of AFP alignments the groups of residues (the group size is the order of symmetry)
	 * that align together and are consistent (equivalent for each subunit). As a result a modified AFPChain with 
	 * <order of symmetry> consistent groups is produced.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] array of the protein.
	 * OUTPUT: an AFP alignment with subunit groups consistent between each other.
	 * ALGORITHM: from the list of AFP alignments an AFP graph that relates the aligned residues is build and for every
	 *            residue in the graph it checks whether there is a cycle of size <order> by performing a DFS of order levels 
	 *            and checking if is possible to return to the first vertex with <order> edges.
	 * RUNNING TIME: order of complexity depending on the order of symmetry and the number of residues (DFS of <order> depth).
	 */
	public static AFPChain refinedAFP(List<int[][][]> allAlignments, Atom[] ca1) throws StructureException {
		
		//Create the alignment graph and initialize a variable to store the groups
		List<List<Integer>> graph = buildAFPgraph(allAlignments, ca1);
		List<List<Integer>> groups = new ArrayList<List<Integer>>();
		
		//Initialize a variable to store the residues already in a group (do not take them twice)
		List<Integer> alreadySeen = new ArrayList<Integer>();
		int order = allAlignments.size()+1;
		
		//Loop through all the residues (vertices) in the graph
		for (int i=0; i<graph.size(); i++){
			if (!alreadySeen.contains(i)){
				
				//Initialize the variables for the DFS of this residue iteration
				Stack<Integer> path = new Stack<Integer>(); //stack that stores the current path nodes
				Stack<List<Integer>> stack = new Stack<List<Integer>>(); //stack that stores the nodes to be visited next: [vertex,level]
				List<Integer> source = new ArrayList<Integer>(); //source information: level 0
				source.add(i);
				source.add(0);
				stack.push(source);
				boolean foundGroup = false; //Do not loop at any other node if you already found a group of connected nodes
				
				while (!stack.isEmpty() && !foundGroup){
					List<Integer> vertex = stack.pop();
					
					//If the vertex level is lower than the past level remove the last element of the path
					if (vertex.get(1)<=path.size() && path.size()>0){
						path.pop();
					}
					
					//If the vertex has level lower than the order consider its neighbors
					if (vertex.get(1)<order){
						//First add the node to the path
						path.push(vertex.get(0));
						
						for (int k=0; k<graph.get(vertex.get(0)).size(); k++){
							//Only add to the stack the nodes not included in the current path
							if (!path.contains(graph.get(vertex.get(0)).get(k)) && !alreadySeen.contains(graph.get(vertex.get(0)).get(k))){
								List<Integer> node = new ArrayList<Integer>();
								node.add(graph.get(vertex.get(0)).get(k));
								node.add(vertex.get(1)+1);
								stack.push(node);
								//If the level=order and the node is equal to the source store the group of residues in the path
								if (node.get(0)==i && node.get(1)==order){
									List<Integer> group = new ArrayList<Integer>();
									path.push(node.get(0));
									int n = path.size();
									System.out.println("Path size: "+n);
									for (int x=0; x<n; x++){
										int p = path.pop();
										group.add(p);
										alreadySeen.add(p);
									}
									group.sort(null);
									for (int e:group){
										System.out.println(e);
									}
									groups.add(group);
									foundGroup = true;
									path.clear();
									break;
								}
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
			for (int k=0; k<group.size(); k++){
				optAln.get(k).get(0).add(group.get(k));
				optAln.get(k).get(1).add(group.get((k+1)%group.size()));
			}
		}
		
		return createOptAln(optAln, ca1);
	}
	
}