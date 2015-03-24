package org.biojava.nbio.structure.align.symm.subunit;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
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
	 * NOT USED ANYMORE.
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
	 * NOT USED ANYMORE, dealing with arrays of AFPChains.
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
		ca2block = Arrays.copyOfRange(ca2, afpChain.getOptAln()[0][1][0], afpChain.getOptAln()[0][1][afpChain.getOptAln()[0][1].length-1]+1);
		
		/*//Try the method in AlignmentTools
		int[] aligned1 = afpChain.getOptAln()[0][0];
		int[] aligned2 = afpChain.getOptAln()[0][1];
		AFPChain displayAFP = AlignmentTools.createAFPChain(ca1block, ca2block, aligned1, aligned2);*/
		
		//Modify the optimal alignment to include only one subunit (block)
		int[][][] optAln = new int[1][2][afpChain.getOptLen()[0]];
		int[][] block = afpChain.getOptAln()[0];
		//Normalize the residues of the second subunit, to be in the range of ca2block
		int start = block[1][0];
		for (int i=0; i<block[1].length; i++){
			block[1][i] -= start;
		}
		optAln[0] = block;
		int[] optLens = new int[1];
		optLens[0]=optAln[0][0].length;
		
		//Modify the AFP chain to adapt the new optimal alignment of two subunits.
		AFPChain displayAFP = new AFPChain();
		try {
			displayAFP = AlignmentTools.replaceOptAln(optAln, afpChain, ca1block, ca2block);
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
			jmolPanel.evalString("select *; backbone off; cartoon on; hide ligand; center;");
			
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
	 * NOT USED ANYMORE, the optAln is an triple array.
	 * 
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
		copyAFP.setBlockGap(AlignmentTools.calculateBlockGap(optAlgn));
		
		//Recalculate properties: superposition, tm-score, etc
		Atom[] ca1clone = StructureTools.cloneCAArray(ca1); // don't modify ca1 positions
		AlignmentTools.updateSuperposition(copyAFP, ca1, ca1clone);
		
		//It re-does the sequence alignment strings from the OptAlgn information only
		copyAFP.setAlnsymb(null);
		AFPAlignmentDisplay.getAlign(copyAFP, ca1, ca1clone);
		
		return copyAFP;
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
	
}