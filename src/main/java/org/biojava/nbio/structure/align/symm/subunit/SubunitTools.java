package org.biojava.nbio.structure.align.symm.subunit;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Stack;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.RotationAxis;

/**
 * Methods to process AFPChain multiple alignments and analyze the symmetrical subunits generated.
 * 
 * @author Aleix Lafita
 * 
 * Last modified: 10.03.2015
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
	 * OUTPUT: a list containing the cutting intervals (size: 2*order), lower and upper bounds of regions.
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
		Collections.sort(intervals);
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
		Collections.sort(selectedInter);
		
		/*//Print the intervals
		System.out.println("Number of interval residues: "+selectedInter.size());
		for (int intr:selectedInter){
			System.out.println(intr);
		}*/
		
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
		Collections.sort(intersectionResidues);
		return intersectionResidues;
	}
	
	/**
	 * Method that processes an AFP alignment and returns a triple list of aligned residues.
	 * 
	 * INPUT: an list AFP alignment.
	 * OUTPUT: a triple list: first index = alignment number, second = chain number, third = position in the alignment.
	 */
	public static List<List<List<Integer>>> processMultipleAFP(List<AFPChain> allAlignments){
		
		//Initialize the triple list to be returned
		List<List<List<Integer>>> alignments = new ArrayList<List<List<Integer>>>();
		
		//Loop through all the alignments and extract a double list of residue pairs
		for (int k=0; k<allAlignments.size(); k++){
			
			//This list contains a list for each chain of the alignments with its residues
			List<List<Integer>> alignment = new ArrayList<List<Integer>>();
			
			for (int x=0; x<2; x++){
				//Get all the residues in one chain of the alignment and add them to the alignment list
				ArrayList<Integer> residues = new ArrayList<Integer>();
				for (int y=0; y<allAlignments.get(k).getOptAln().length; y++){
					
					//Get the total number of residues in the group
					int n = allAlignments.get(k).getOptAln()[y][x].length;
					//Append every residue to a list of residues
					for (int j=0; j<n; j++){
						residues.add(allAlignments.get(k).getOptAln()[y][x][j]);
					}
				}
				alignment.add(residues);
			}
			alignments.add(alignment);
		}
		return alignments;
	}
	
	/**
	 * Method that displays two superimposed subunits in jmol.
	 * 
	 * INPUT: an AFP alignment and the protein data.
	 * OUTPUT: a jmol panel with only one subunit superimposed.
	 */
	public static void displaySuperimposedSubunits(AFPChain afpChain, String name, Atom[] ca1, Atom[] ca2){
		
		//Create the atom arrays corresponding to the first and second subunits only
		Atom[] ca1block = new Atom[afpChain.getOptLen()[0]];
		Atom[] ca2block = new Atom[afpChain.getOptLen()[0]];
		ca1block = Arrays.copyOfRange(ca1, 0, afpChain.getOptAln()[0][0][afpChain.getOptAln()[0][0].length-1]+1);
		ca2block = Arrays.copyOfRange(ca2, 0, afpChain.getOptAln()[0][1][afpChain.getOptAln()[0][1].length-1]+1);
		
		//Modify the optimal alignment to include only one subunit (block)
		AFPChain displayAFP = (AFPChain) afpChain.clone();
		int[][][] optAln = new int[1][2][displayAFP.getOptLen()[0]];
		int[][] block1 = displayAFP.getOptAln()[0];
		optAln[0] = block1;
		int[] optLens = new int[1];
		optLens[0]=optAln[0][0].length;
		
		//Modify the AFP chain to adapt the new optimal alignment of two subunits.
		try {
			displayAFP = SubunitTools.createOptAln(optAln, displayAFP, ca1block, ca2block);
		} catch (StructureException e1) {
			e1.printStackTrace();
		}
		
		//Another array to display is created only with the residues of the second subunit, because all (first and second, are needed to superimpose, but only the second is relevant in the alignment)
		//DOES NOT WORK, because the second subunit is not colored
		//Atom[] ca2blockDisplay = Arrays.copyOfRange(ca2, afpChain.getOptAln()[0][1][0], afpChain.getOptAln()[0][1][afpChain.getOptAln()[0][1].length-1]+1);

		//Set the name of the protein
		displayAFP.setName1(name+" su1");
		displayAFP.setName2(name+" su2");
		
		try {
			
			//Display the AFP alignment of the subunits
			StructureAlignmentJmol jmolPanel;
			jmolPanel = StructureAlignmentDisplay.display(displayAFP, ca1block, ca2block);
			jmolPanel.evalString("select *; backbone off; cartoon on; model 1; hide ligand; center;");
			
			/*	
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(displayAFP);
			jmolPanel.evalString(axis.getJmolScript(ca1));*/
			
		} catch (StructureException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Method that displays the structural alignment with each subunit colored differently.
	 * 
	 * INPUT: an AFP alignment and the protein data.
	 * OUTPUT: displays the protein alignment with jmol coloring the subunit residues.
	 */
	public static void displayColorSubunits(AFPChain afpChain, String name, Atom[] ca1, Atom[] ca2){
		
		//Set the name of the protein
		afpChain.setName1(name);
		afpChain.setName2(name);
		
		try {
			
			//Display the AFP alignment of the subunits
			StructureAlignmentJmol jmolPanel;
			jmolPanel = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
				
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(afpChain);
			jmolPanel.evalString(axis.getJmolScript(ca1));
			
		} catch (StructureException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * NOT USED ANYMORE: new method refinedAFP uses only a triple list as an input, instead of the AFP alignments.
	 * 
	 * It returns an AFPChain with a segmented optimal alignment, which means that it has <order of symmetry> blocks of
	 * aligned subunits.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] arrays.
	 * OUTPUT: the optimal AFPChain alignment divided into the subunits.
	 */
	public static AFPChain createOptAlgn(List<AFPChain> allAlignments, Atom[] ca1) throws StructureException {
		
		//Extract the subunits from the set of alignments
		List<List<Integer>> subunits = extractSubunits(ca1, allAlignments);
		int order = subunits.size();
		
		AFPChain optimalAFP = allAlignments.get(0);
				
		//Extract the old alignment
		int[][][] oldAlgn = optimalAFP.getOptAln();
		
		//Create a triple list with the optimal alignment divided into the subunits in <order> blocs
		List<List<List<Integer>>> newAlgn = new ArrayList<List<List<Integer>>>();
		
		//Loop through all the optimal alignment and subunits and select the pairs of aligned residues for each subunit
		for (int i=0; i<order; i++){
			List<Integer> chain1 = new ArrayList<Integer>();
			List<Integer> chain2 = new ArrayList<Integer>();
			for (int k=0; k<oldAlgn.length; k++){
				for (int j=0; j<oldAlgn[k][0].length; j++){
					if (subunits.get(i).contains(oldAlgn[k][0][j])){
						chain1.add(oldAlgn[k][0][j]);
						chain2.add(oldAlgn[k][1][j]);
					}
				}
			}
			List<List<Integer>> subunit = new ArrayList<List<Integer>>();
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
		AFPChain refinedAFP = (AFPChain) allAlignments.get(0).clone();
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
		
		//It re-does the sequence alignment strings from the OptAlgn information only
		AFPAlignmentDisplay.getAlign(refinedAFP, ca1, ca1);
				
		return refinedAFP;
	}
	
	/**
	 * It returns an AFPChain with a segmented optimal alignment, which means that it has <order of symmetry> blocks of
	 * aligned subunits. The original AFPChain is not modified.
	 * 
	 * INPUT: the optimal alignment in a triple list (same format as optAln of AFPChain) and the Atom[] arrays.
	 * OUTPUT: the optimal AFPChain alignment object divided into the subunits.
	 */
	public static AFPChain createOptAln(int[][][] newAlgn, AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException {
		
		//The order is the number of groups in the newAlgn
		int order = newAlgn.length;
		
		//Calculate the alignment length from all the subunits lengths
		int[] optLens = new int[order];
		for(int s=0;s<order;s++) {
			optLens[s] = newAlgn[s][0].length;
		}
		int optLength = 0;
		for(int s=0;s<order;s++) {
			optLength += optLens[s];
		}
		
		//Print the sizes to check correctness
		System.out.println("Number of subunits: "+newAlgn.length);
		System.out.println("Subunit length: "+newAlgn[0][0].length);
		
		//Create a copy of the original AFPChain and set everything needed for the structure update
		AFPChain copyAFP = (AFPChain) afpChain.clone();
		
		//Set the new parameters of the optimal alignment
		copyAFP.setOptLength(optLength);
		copyAFP.setOptLen(optLens);
		copyAFP.setOptAln(newAlgn);
		
		//Set the block information of the new alignment
		copyAFP.setBlockNum(order);
		copyAFP.setBlockSize(optLens);
		copyAFP.setBlockResList(newAlgn);
		copyAFP.setBlockResSize(optLens);
		copyAFP.setBlockGap(createBlockGap(newAlgn));
		
		//Recalculate properties: superposition, tm-score, etc
		Atom[] ca2clone = StructureTools.cloneCAArray(ca2); // don't modify ca1 positions
		AlignmentTools.updateSuperposition(copyAFP, ca1, ca2clone);
		
		//It re-does the sequence alignment strings from the OptAlgn information only
		copyAFP.setAlnsymb(null);
		AFPAlignmentDisplay.getAlign(copyAFP, ca1, ca2clone);
		
		return copyAFP;
	}
	
	/**
	 * It returns an AFPChain with a segmented optimal alignment, which means that it has <order of symmetry> blocks of
	 * aligned subunits. The original AFPChain is not modified.
	 * 
	 * INPUT: the optimal alignment in a triple list (same format as optAln of AFPChain) and the Atom[] arrays.
	 * OUTPUT: the optimal AFPChain alignment object divided into the subunits.
	 */
	public static AFPChain createOptAln(List<List<List<Integer>>> newAlgn, AFPChain afpChain, Atom[] ca1) throws StructureException {
		
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
		
		//Print the sizes to check correctness
		System.out.println("Number of subunits: "+optAlgn.length);
		System.out.println("Subunit length: "+optAlgn[0][0].length);
		
		//Create a copy of the original AFPChain and set everything needed for the structure update
		AFPChain copyAFP = (AFPChain) afpChain.clone();
		
		//Set the new parameters of the optimal alignment
		copyAFP.setOptLength(optLength);
		copyAFP.setOptLen(optLens);
		copyAFP.setOptAln(optAlgn);
		
		//Set the block information of the new alignment
		copyAFP.setBlockNum(order);
		copyAFP.setBlockSize(optLens);
		copyAFP.setBlockResList(optAlgn);
		copyAFP.setBlockResSize(optLens);
		copyAFP.setBlockGap(createBlockGap(optAlgn));
		
		//Recalculate properties: superposition, tm-score, etc
		Atom[] ca1clone = StructureTools.cloneCAArray(ca1); // don't modify ca1 positions
		AlignmentTools.updateSuperposition(copyAFP, ca1, ca1clone);
		
		//It re-does the sequence alignment strings from the OptAlgn information only
		copyAFP.setAlnsymb(null);
		AFPAlignmentDisplay.getAlign(copyAFP, ca1, ca1clone);
		
		return copyAFP;
	}
	
	/**
	 * Calculates a graph in the format of adjacency list from the set of alignments, where each vertex is a 
	 * residue and each edge means the connection between the two residues in one of the alignments.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] array.
	 * OUTPUT: the alignment graph (describing relations between residues). List dimensions: AdjList[vertices][edges]
	 */
	public static List<List<Integer>> buildAFPgraph(List<AFPChain> allAlignments, Atom[] ca1) {
		
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
				/*//Make the graph undirected (optional feature)
				if (!adjList.get(edge).contains(vertex)){
					adjList.get(edge).add(vertex);
				}*/
			}
		}
		//Print graph information
		System.out.println("GRAPH INFORMATION:");
		System.out.println(" - Vertices: "+adjList.size());
		int edges = 0;
		for (List<Integer> v:adjList){
			edges+=v.size();
		}
		System.out.println(" - Edges: "+edges);
		//System.out.println(" - Consistent connected residues: "+count+" out of "+adjList.size());
		
		//Store the graph as a csv file to analyze it graphically
		//generateCsvFile(adjList, "unknown_dir");
		
		return adjList;
	}
	
	/**
	 * Calculates a weighted graph in the format of a matrix from the set of alignments, where each vertex is a 
	 * residue and each edge means the connection between the two residues in one of the alignments. The weight 
	 * of the edge is the distance between both residues in the protein.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] array.
	 * OUTPUT: the alignment graph (describing relations between residues). List dimensions: AdjList[vertices][edges]
	 */
	public static double[][][] buildWeightedAFPgraph(List<AFPChain> allAlignments, Atom[] ca1) {
		
		//Initialize the matrix that stores the graph and fill it with 0 values
		double[][][] graph = new double[ca1.length][][];
		for (int i=0; i<ca1.length; i++){
			double[][] row = new double[ca1.length][];
			for (int j=0; j<ca1.length; j++){
				double[] entry = {0.0,10.0};
				row[j] = entry;
			}
			graph[i] = row;
		}
		
		//Convert the multiple alignments into a triple list to handle the relations
		List<List<List<Integer>>> alignments = processMultipleAFP(allAlignments);
		
		for (int k=0; k< alignments.size(); k++){
			for (int i=0; i<alignments.get(k).get(0).size(); i++){
				
				//The vertex is the residue in the first chain and the edge the one in the second chain
				int vertex = alignments.get(k).get(0).get(i);
				int edge = alignments.get(k).get(1).get(i);
				if (graph[vertex][edge][0]==0.0 && graph[edge][vertex][0]==0.0){
					//Distance between two residues in the protein
					graph[vertex][edge][0] = Calc.getDistance(ca1[vertex],ca1[edge]);
					//Distance between two residues in the superimposed alignment
					graph[vertex][edge][1] = allAlignments.get(k).getDistanceMatrix().getArray()[vertex][edge];
				}
			}
		}
		
		//Store the graph as a csv file to analyze it graphically
		csvWeightedGraph(graph, "unknown_align");
		
		return graph;
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
	 *            and checking if is possible to return to the first vertex with <order> edges. Because this selection is 
	 *            very restrictive and does not select some of the symmetric residues, a last step selecting the residues 
	 *            with <order> neighbors and its neighbors if they are also consistent can be done.
	 * RUNNING TIME: order of complexity depending on the order of symmetry and the number of residues (DFS of <order> depth).
	 *               The running time is much faster than the multiple alignment so it is not a bottleneck.
	 */
	public static AFPChain refinedAFPold(List<AFPChain> allAlignments, Atom[] ca1) throws StructureException {
		
		//Create the alignment graph and initialize a variable to store the groups
		List<List<Integer>> graph = buildAFPgraph(allAlignments, ca1);
		List<List<Integer>> groups = new ArrayList<List<Integer>>();
		List<Integer> intervals = calculateIntervals(ca1, allAlignments);
		//Add to the intervals list the first and last element of the protein to consider also cycles in this region (permissive intervals)
		intervals.add(0, 0);
		intervals.add(ca1.length-1);
		
		//Initialize a variable to store the residues already in a group (do not take them twice)
		List<Integer> alreadySeen = new ArrayList<Integer>();
		int order = allAlignments.size()+1;
		
		List<Integer> lastGroup = new ArrayList<Integer>();
		for (int i=0; i<order; i++){
			lastGroup.add(0);
		}
		
		//Loop through all the residues (vertices) in the graph
		for (int i=0; i<graph.size(); i++){
			if (!alreadySeen.contains(i)){
				//System.out.println("Cycle for residue "+i);
				
				//Initialize the variables for the DFS of this residue iteration
				Stack<Integer> path = new Stack<Integer>(); //stack that stores the current path nodes
				Stack<List<Integer>> stack = new Stack<List<Integer>>(); //stack that stores the nodes to be visited next: [vertex,level]
				List<Integer> source = new ArrayList<Integer>(); //source information: level 0
				source.add(i);
				source.add(0);
				stack.push(source);
				
				boolean foundGroup = false; //Do not loop at any other node if you already found a group of connected nodes
				
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
					if (vertex.get(1)<=path.size()-1){
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
							if (!path.contains(node.get(0)) && !alreadySeen.contains(node.get(0)) && node.get(1)<order){
								stack.push(node);
							}
							//If the level=order and the node is equal to the source store the group of residues in the path
							else if (node.get(0)==i && node.get(1)==order){

								List<Integer> group = new ArrayList<Integer>();
								int n = path.size();
								//System.out.println("Path size: "+n);
								for (int x=0; x<n; x++){
									int p = path.get(x);
									group.add(p);
								}
								Collections.sort(group);
								boolean correct = true;
								
								//Check that all the residues are greater than the last group found, otherwise inconsistent
								for (int e=0; e<group.size()-1; e++){
									if (!(group.get(e)>lastGroup.get(e)) && lastGroup.get(e+1)!=0){
										//System.out.println("Inconsistent group: not sequential order...");
										//System.out.println("Last group "+lastGroup.get(e)+", Current group "+group.get(e));
										correct = false;
										break;
									}
								}
								//Check that all the residues are lower than the next residue of the last group found
								for (int e=0; e<group.size()-1; e++){
									if (group.get(e)>lastGroup.get(e+1) && lastGroup.get(e+1)!=0 && correct){
										//System.out.println("Inconsistent group: not consistent with the last group...");
										//System.out.println("Last group "+lastGroup.get(e+1)+", Current group "+group.get(e));
										correct = false;
										break;
									}
								}
								//Also check that the residues are inside its subunit intervals boundaries (restrictive intervals)
								for (int e=0; e<group.size(); e++){
									if (!((group.get(e)>=intervals.get(2*e+1)) && (group.get(e)<=intervals.get(2*e+2))) && correct){
										//System.out.println("Inconsistent group: not inside interval boundaries...");
										//System.out.println("Residue "+group.get(e)+" not in the interval ["+intervals.get(2*e+1)+", "+intervals.get(2*e+2)+"]");
										correct = false;
										break;
									}
								}
								if (correct){
									/*System.out.println("Group size: "+group.size());
									for (int e=0; e<group.size(); e++){
										System.out.println(group.get(e));
									}*/
									groups.add(group);
									//Add the vertices of the path to alreadySeen to avoid including them again later groups
									for (int p:group){
										alreadySeen.add(p);
									}
									foundGroup = true;
									path.clear();
									//System.out.println("Clearing path...");
									lastGroup = group;
									break;
								}
							}
						}
					}
				} //end of DFS
			}
		} //end of all the residue analysis
		
		//Last step to select those residues that do not form a cycle but that are connected in a consistent manner.
		boolean last_step = false; //set to true if this last step is required
		if (last_step){
			for (int i=0; i<graph.size(); i++){
				if (graph.get(i).size()==order-1 && !alreadySeen.contains(i)){
					boolean correct = true;
					List<Integer> group = new ArrayList<Integer>();
					//Check that all the residues in the group have not been already seen
					group.add(i);
					for (int j=0; j<graph.get(i).size(); j++){
						group.add(graph.get(i).get(j));
						if (alreadySeen.contains(group.get(j))){
							correct = false;
						}
					}
					Collections.sort(group);
					//System.out.println("Considering a group...");
					/*for (int e=0; e<group.size(); e++){
						System.out.println(group.get(e));
					}*/
					
					//Check that the residues are inside its subunit intervals boundaries (Needed to avoid discontinuous regions) (restrictive intervals)
					for (int e=0; e<group.size(); e++){
						if (!((group.get(e)>=intervals.get(2*e+1)) && (group.get(e)<=intervals.get(2*e+2))) && correct){
							//System.out.println("Inconsistent group: not inside interval boundaries...");
							//System.out.println("Residue "+group.get(e)+" not in the interval ["+intervals.get(2*e)+", "+intervals.get(2*e+3)+"]");
							correct = false;
							break;
						}
					}
					//Find the insertion index of the group in the groups list and check that all other residues are consistent
					int index = 0;
					for (int k=0; k<groups.size(); k++){
						//special case for the first group, compare only with the first
						if (k==0 && groups.get(k).get(0)>group.get(0)){
							index = k;
							for (int d=0; d<order; d++){
								if (groups.get(k).get(d)<group.get(d)){
									correct = false;
									//System.out.println("Inconsistent group: not accordance with neighbors...");
									//System.out.println("Residue "+group.get(d)+" not in the neighbors interval [0, "+groups.get(k).get(d)+"]");
									break;
								}
							}
							break;
						}
						//case for all other values of k, compare previous and next
						else if (groups.get(k).get(0)>group.get(0)){
							index = k;
							for (int d=0; d<order; d++){
								if (groups.get(k-1).get(d)>group.get(d) || groups.get(k).get(d)<group.get(d)){
									correct = false;
									//System.out.println("Inconsistent group: not accordance with neighbors...");
									//System.out.println("Residue "+group.get(d)+" not in the neighbors interval ["+groups.get(k-1).get(d)+", "+groups.get(k).get(d)+"]");
									break;
								}
							}
							break;
						}
						if (!correct){
							break;
						}
					}
					//If the conditions are fulfilled insert the group into the right position and mark the residues as seen
					if (correct){
						if (index==0){
							groups.add(group);
						}
						else{
							groups.add(index, group);
						}
						//System.out.println("Extra group added: ");
						for (int e:group){
							alreadySeen.add(e);
							//System.out.println(e);
						}
					}
				}
			}
		}
		
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
				optAln.get(k).get(1).add(group.get((k+1)%order));
				optAln.get(k).get(0).add(group.get(k));
			}
		}
		return createOptAln(optAln, allAlignments.get(order-2), ca1);
	}
	
	/**
	 * Method that calculates the number of gaps in each subunit block of an optimal AFP alignment.
	 * 
	 * INPUT: an optimal alignment in the format int[][][].
	 * OUTPUT: an int[] array of <order> length containing the gaps in each block as int[block].
	 */
	public static int[] createBlockGap(int[][][] optAln){
		
		//Initialize the array to be returned
		int [] blockGap = new int[optAln.length];
		
		//Loop for every block and look in both chains for non-contiguous residues.
		for (int i=0; i<optAln.length; i++){
			int gaps = 0; //the number of gaps in that block
			int last1 = 0; //the last residue position in chain 1
			int last2 = 0; //the last residue position in chain 2
			//Loop for every position in the block
			for (int j=0; j<optAln[i][0].length; j++){
				//If the first position is evaluated initialize the last positions
				if (j==0){
					last1 = optAln[i][0][j];
					last2 = optAln[i][1][j];
				}
				else{
					//If one of the positions or both are not contiguous increment the number of gaps
					if (optAln[i][0][j] > last1+1 || optAln[i][1][j] > last2+1){
						gaps++;
						last1 = optAln[i][0][j];
						last2 = optAln[i][1][j];
					}
					//Otherwise just set the last position to the current one
					else{
						last1 = optAln[i][0][j];
						last2 = optAln[i][1][j];
					}
				}
			}
			blockGap[i] = gaps;
		}
		return blockGap;
	}
	
	
	/**
	 * Saves a graph into a csv file in the format of tuples (vertex,edge) for every edge in the graph.
	 */
	 public static void csvGraph(List<List<Integer>> graph, String name){
		 
		String sFileName = "/home/scratch/graphs/"+name+".csv";
		
		try
		{
		    FileWriter writer = new FileWriter(sFileName);
		    writer.append("Vertex,Edge\n");
		    for (Integer i=0; i<graph.size(); i++){
		    	for (int j=0; j<graph.get(i).size(); j++){
		    		
		    		writer.append(i.toString());
		    		writer.append(',');
		    		writer.append(graph.get(i).get(j).toString());
		    		writer.append('\n');
		    	}
		    }
		    
		    writer.flush();
		    writer.close();
		}
		catch(IOException e)
		{
		     e.printStackTrace();
		}
	}
	 
	/**
	 * Saves a graph into a csv file in the format of tuples (vertex,edge) for every node in the graph that is part of a subunit.
	 */
	 public static void csvGraphSubunits(List<List<Integer>> graph, String name, List<Integer> alreadySeen){
		 
		String sFileName = "/home/scratch/graphs/"+name+"_subunit.csv";
		
		try
		{
		    FileWriter writer = new FileWriter(sFileName);
		    writer.append("Vertex,Edge\n");	    
		    for (int i=0; i<graph.size(); i++){
		    	if (alreadySeen.contains(i)){
			    	for (int j=0; j<graph.get(i).size(); j++){
			    		if (alreadySeen.contains(graph.get(i).get(j))){
			    		writer.append(i+",");
			    		writer.append(graph.get(i).get(j)+"\n");
			    		}
			    	}
		    	}
		    }
		    
		    writer.flush();
		    writer.close();
		}
		catch(IOException e)
		{
		     e.printStackTrace();
		}
	}
	
	 /**
	 * Saves a graph into a csv file in the format of tuples (vertex,edge,weight) for every edge in the graph.
	 */
	public static void csvWeightedGraph(double[][][] graph, String name){
		 
		String sFileName = "/home/scratch/graphs/"+name+".csv";
		
		try
		{
		    FileWriter writer = new FileWriter(sFileName);
		    writer.append("Vertex,Edge,Distance,RMSD\n");	    
		    for (int i=0; i<graph.length; i++){
		    	for (int j=0; j<graph[i].length; j++){
		    		
		    		if (graph[i][j][0]!=0.0){
			    		writer.append(i+",");
			    		writer.append(j+",");
			    		writer.append(graph[i][j][0]+",");
			    		writer.append(graph[i][j][1]+"\n");
		    		}
		    	}
		    }
		    writer.flush();
		    writer.close();
		}
		catch(IOException e)
		{
		     e.printStackTrace();
		}
	}
	 
	 /**
	 * Calculates from a set of AFP alignments the groups of residues (the group size is the order of symmetry)
	 * that align together and are consistent, equivalent for each subunit). As a result a modified AFPChain with 
	 * <order of symmetry> consistent groups is produced.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] array of the protein.
	 * OUTPUT: an AFP alignment with subunit groups consistent between each other.
	 * ALGORITHM: TBD
	 * 
	 * RUNNING TIME: the order of complexity is polynomial in length with the exponent being the order of symmetry.
	 */
	public static AFPChain refinedAFP(List<AFPChain> allAlignments, Atom[] ca1) throws StructureException {
		
		//Create the alignment graph and initialize a variable to store the groups
		double[][][] weights = buildWeightedAFPgraph(allAlignments, ca1);
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
		int order = allAlignments.size()+1;
		
		//Loop through all the residues (vertices) in the graph
		for (int i=0; i<graph.size(); i++){
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
							
							//Check that the residues have consistent interresidue distance between the group. The total distance should be more or less the same
							double[] distances = new double[order];
							double maxDist = 0.0;
							boolean consistent = true;
							
							for (int g=0; g<order; g++){
								for (int h=0; h<order; h++){
									distances[g] += Calc.getDistance(ca1[group.get(g)],ca1[group.get(h)]);
									int residueDist = Math.abs(group.get(g)-group.get(h));
									//If the distance between the two residues in number is lower they are too close and not consistent.
									if (residueDist<(aligned_res/(order+2)) && h!=g){
										consistent = false;
										//System.out.println("Not consistent: difference of "+residueDist+" with maximum "+aligned_res/(order+1)+"...");
										break;
									}
								}
								if (maxDist<distances[g])
									maxDist=distances[g];
							}
							
							for (int g=0; g<order; g++){
								for (int h=0; h<order; h++){
									//If their difference in distance is higher than 15% they are not consistent
									if (Math.abs(distances[g]-distances[h])>maxDist/order){
										consistent = false;
										//System.out.println("Not consistent: difference of "+Math.abs(distances[g]-distances[h])+" with maximum "+maxDist+"...");
										break;
									}
								}
							}
							if (!consistent) continue;
							
							//If any residue of the group is in another group already seen mark it as inconsistent
							for (int d=0; d<order; d++){
								if (alreadySeen.contains(group.get(d))){
									consistent = false;
								}
							}
							
							boolean added = false;
							//Check that the group is consistent with the groups already found
							//Find the insertion index of the group in the groups list and replace another group if necessary
							int index = -1;
							for (int j=0; j<groups.size(); j++){
								//special case for the first group, compare only with the first
								if (j==0 && groups.get(j).get(0)>=group.get(0)){
									index = j;
									for (int d=0; d<order; d++){
										if (groups.get(j).get(d)<=group.get(d)){
											//Compare both groups and decide which is more meaningful based on their RMSD (the lower wins)
											if (groups.get(j).get(d)==group.get(d)){
												double current_rmsd = 0.0;
												double new_rmsd = 0.0;
												for (int g=0; g<order; g++){
													for (int h=0; h<order; h++){
														current_rmsd += weights[groups.get(j).get(g)][groups.get(j).get(h)][1];
														new_rmsd += weights[group.get(g)][group.get(h)][1];
													}
												}
												if (new_rmsd<current_rmsd){
													groups.remove(j);
													groups.add(j, group);
													added = true;
													foundGroup = true;
													//System.out.println("Replaced the group of residue "+j);
												}
												else{
													consistent=false;
												}
											}
											break;
										}
									}
									break;
								}
								//special case for the last group, compare only with the last
								else if (j==group.size()-1 && groups.get(j).get(0)<=group.get(0)){
									index = -1;
									for (int d=0; d<order; d++){
										if (groups.get(j).get(d)>=group.get(d)){
											//Compare both groups and decide which is more meaningful based on their RMSD (the lower wins)
											if (groups.get(j).get(d)==group.get(d)){
												double current_rmsd = 0.0;
												double new_rmsd = 0.0;
												for (int g=0; g<order; g++){
													for (int h=0; h<order; h++){
														current_rmsd += weights[groups.get(j).get(g)][groups.get(j).get(h)][1];
														new_rmsd += weights[group.get(g)][group.get(h)][1];
													}
												}
												if (new_rmsd<current_rmsd){
													groups.remove(j);
													groups.add(j, group);
													added = true;
													foundGroup = true;
													//System.out.println("Replaced the group of residue "+j);
												}
												else{
													consistent=false;
													//System.out.println("Not replaced the group of residue "+j);
												}
											}
											break;
										}
									}
									break;
								}
								//case for all other values of k, compare previous and next
								else if (groups.get(j).get(0)>group.get(0)){
									index = j;
									for (int d=0; d<order; d++){
										if (groups.get(j).get(d)<=group.get(d) || groups.get(j-1).get(d)>=group.get(d)){
											//Compare both groups and decide which is more meaningful based on their RMSD (the lower wins)
											if (groups.get(j).get(d)==group.get(d)){
												double current_rmsd = 0.0;
												double new_rmsd = 0.0;
												for (int g=0; g<order; g++){
													for (int h=0; h<order; h++){
														current_rmsd += weights[groups.get(j).get(g)][groups.get(j).get(h)][1];
														new_rmsd += weights[group.get(g)][group.get(h)][1];
													}
												}
												if (new_rmsd<current_rmsd){
													groups.remove(j);
													groups.add(j, group);
													added = true;
													foundGroup = true;
													//System.out.println("Replaced the group of residue "+j);
												}
												else{
													consistent=false;
													//System.out.println("Not replaced the group of residue "+j);
												}
											}
											else if (groups.get(j-1).get(d)==group.get(d)){
												double current_rmsd = 0.0;
												double new_rmsd = 0.0;
												for (int g=0; g<order; g++){
													for (int h=0; h<order; h++){
														current_rmsd += weights[groups.get(j-1).get(g)][groups.get(j-1).get(h)][1];
														new_rmsd += weights[group.get(g)][group.get(h)][1];
													}
												}
												if (new_rmsd<current_rmsd){
													groups.remove(j-1);
													groups.add(j-1, group);
													added = true;
													foundGroup = true;
													//System.out.println("Replaced the group of residue "+(j-1));
													//System.out.println("Not replaced the group of residue "+j);
												}
												else{
													consistent=false;
												}
											}
											break;
										}
									}
									break;
								}
							}
							if (!added && consistent){
								//If the conditions are fulfilled insert the group into the right position
								if (index==-1){
									groups.add(group);
								}
								else{
									groups.add(index, group);
								}
								System.out.println("Group added, size: "+group.size());
								for (int e:group){
									alreadySeen.add(e);
									System.out.println(e);
								}
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
		
		//Save the graph with only the subunit residues, to see the filtering made by the algorithm
		csvGraphSubunits(graph, "unknown", alreadySeen);
		
		return createOptAln(optAln, allAlignments.get(order-2), ca1);
	}
}