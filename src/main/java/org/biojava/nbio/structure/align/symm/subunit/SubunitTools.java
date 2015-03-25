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
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;

/**
 * Methods to process AFPChain multiple alignments and analyze the symmetrical subunits generated.
 * 
 * @author lafita
 */
public class SubunitTools {
	
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
			
			/*	
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(displayAFP);
			jmolPanel.evalString(axis.getJmolScript(ca1));*/
			
		} catch (StructureException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates a weighted graph in the format of a matrix from the set of alignments, where each vertex is a 
	 * residue and each edge means the connection between the two residues in one of the alignments. The weight 
	 * of the edge is the distance between both residues in the protein.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] array.
	 * OUTPUT: the alignment graph (describing relations between residues). List dimensions: AdjList[vertices][edges]
	 */
	public static double[][][] buildWeightedAFPgraph(AFPChain[] allAlignments, Atom[] ca1) throws IOException {
		
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
		
		for (int k=0; k< allAlignments.length; k++){
			for (int i=0; i<allAlignments[k].getOptAln().length; i++){
				for (int j=0; j<allAlignments[k].getOptAln()[i][0].length; j++){
					//The vertex is the residue in the first chain and the edge the one in the second chain
					int vertex = allAlignments[k].getOptAln()[i][0][j];
					int edge = allAlignments[k].getOptAln()[i][1][j];
					if (graph[vertex][edge][0]==0.0 && graph[edge][vertex][0]==0.0){
						//Distance between two residues in the protein
						graph[vertex][edge][0] = Calc.getDistance(ca1[vertex],ca1[edge]);
						//Distance between two residues in the superimposed alignment
						graph[vertex][edge][1] = allAlignments[k].getDistanceMatrix().getArray()[vertex][edge];
					}
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
	 public static void csvGraph(List<List<Integer>> graph, String name) throws IOException {
		 
		String sFileName = "/home/scratch/graphs/"+name+".csv";
		
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
	 
	/**
	 * Saves a graph into a csv file in the format of tuples (vertex,edge) for every node in the graph that is part of a subunit.
	 */
	 public static void csvGraphSubunits(List<List<Integer>> graph, String name, List<Integer> alreadySeen) throws IOException{
		 
	String sFileName = "/home/scratch/graphs/"+name+"_subunit.csv";

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
	
	 /**
	 * Saves a graph into a csv file in the format of tuples (vertex,edge,weight) for every edge in the graph.
	 */
	public static void csvWeightedGraph(double[][][] graph, String name) throws IOException{
		 
		String sFileName = "/home/scratch/graphs/"+name+".csv";
		
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
	
	 /**
	 * Saves a graph into a csv file in the format of tuples (vertex,edge,weight) for every edge in the graph.
	 * @throws StructureException 
	 */
	public static void analyzeRunningTime(String[] names) throws IOException, StructureException{
		
		//Prepare the file writer
		String sFileName = "/home/scratch/stats/complexity.csv";
	    FileWriter writer = new FileWriter(sFileName);
	    writer.append("Name,Order,Length,TimeMultiple,TimeSingle,TimeNotRefined\n");
		
		for (int i=0; i<names.length; i++){
			
			String name = names[i];
			int order = 0;
			int length = 0;
			
			//Initialize the time measure variables
			long durationMultiple = 0;
			long durationSingle = 0;
			long durationNoRefine = 0;
			
			for (int j=0; j<3; j++){
			
				long startTime = System.nanoTime();
				
				//Set the name of the protein structure to analyze
				System.out.println("Analyzing protein "+name );
				AtomCache cache = new AtomCache();
	
				//Parse atoms of the protein into two DS
				Atom[] ca1 = cache.getAtoms(name);
				Atom[] ca2 = cache.getAtoms(name);
				
				//Initialize a new CeSymm class and its parameters and a new alignment class
				CeSymm ceSymm = new CeSymm();
				AFPChain afpChain = new AFPChain();
	
				if (j==0){
					//Perform the alignment
					afpChain = ceSymm.align(ca1, ca2);
					
					long endTime = System.nanoTime();
					
					durationMultiple = (endTime - startTime);
					length = ca1.length;
					order = afpChain.getBlockNum();
				}
				else if(j==1){
					//Perform the alignment
					CESymmParameters params = new CESymmParameters();
					params.setRefineMethod(RefineMethod.SINGLE);
					
					afpChain = ceSymm.align(ca1, ca2, params);
					
					long endTime = System.nanoTime();
					
					durationSingle = (endTime - startTime);
				}
				else{
					//Perform the alignment
					CESymmParameters params = new CESymmParameters();
					params.setRefineMethod(RefineMethod.NOT_REFINED);
					
					afpChain = ceSymm.align(ca1, ca2, params);
					
					long endTime = System.nanoTime();
					
					durationNoRefine = (endTime - startTime);
				}
			}
			writer.append(name+","+order+","+length+","+durationMultiple+","+durationSingle+","+durationNoRefine+"\n");
		}
		
	    writer.flush();
	    writer.close();
	}
	
	public static void main(String[] args) throws StructureException, IOException{
		
		String[] names = {"2F9H.A", "1SQU.A", "3HDP", "2AFG.A", "4DOU", "1VYM", "1HCE", "1TIE", "4I4Q", "1GEN", "1HXN", "1G61.A","1TL2.A","2JAJ.A", "1U6D", "1JOF.A", "1JTD.B",  "1A12.A", "2I5I.A", "1K3I.A", "1GOT.B", "1TIM.A", "1VZW", "1NSJ"};
		analyzeRunningTime(names);
	}
}