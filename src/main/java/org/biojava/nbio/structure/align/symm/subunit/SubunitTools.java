package org.biojava.nbio.structure.align.symm.subunit;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * NOT PART OF THE CE-Symm CORE. Only contais methods to analyze running time and save alignment graphs.
 * 
 * @author lafita
 */
public class SubunitTools {
	
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
	 * Saves a graph into a csv file in the format of tuples (vertex,edge,weight) for every edge in the graph.
	 */
	 public static void csvWeightedGraph(double[][][] graph, String sFileName) throws IOException{
		
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
	 * @throws OrderDetectionFailedException 
	 */
	public static void analyzeRunningTime(String[] names, int[] orders, String sFileName) throws IOException, StructureException, OrderDetectionFailedException{
		
		//Prepare the file writer
	    FileWriter writer = new FileWriter(sFileName);
	    writer.append("Name,Order,orderMultipe,orderSingle,length,timeMultiple,timeSingle,timeNotRefined,timeMultipleOpt\n");
		
		for (int i=0; i<names.length; i++){
			
			String name = names[i];
			int orderM = 0;
			int orderS = 0;
			int length = 0;
			
			//Initialize the time measure variables
			long durationMultiple = 0;
			long durationSingle = 0;
			long durationNoRefine = 0;
			long durationMultipleOpt = 0;
			
			for (int j=0; j<4; j++){
			
				long startTime = System.nanoTime();
				
				//Set the name of the protein structure to analyze
				System.out.println("Analyzing protein "+name);
				AtomCache cache = new AtomCache();
	
				//Parse atoms of the protein into two DS
				Atom[] ca1 = cache.getAtoms(name);
				Atom[] ca2 = cache.getAtoms(name);
				
				//Initialize a new CeSymm class and its parameters and a new alignment class
				CeSymm ceSymm = new CeSymm();
				AFPChain afpChain = new AFPChain();
	
				if (j==0){
					//Perform the alignment
					CESymmParameters params = new CESymmParameters();
					params.setRefineMethod(RefineMethod.MULTIPLE);
					
					afpChain = ceSymm.align(ca1, ca2);
					
					orderM = afpChain.getBlockNum();
					
					long endTime = System.nanoTime();
					
					durationMultiple = (endTime - startTime);
					length = ca1.length;
				}
				else if(j==1){
					//Perform the alignment
					CESymmParameters params = new CESymmParameters();
					params.setRefineMethod(RefineMethod.SINGLE);
					
					afpChain = ceSymm.align(ca1, ca2, params);
					SequenceFunctionOrderDetector orderDetector= new SequenceFunctionOrderDetector();
					
					orderS = orderDetector.calculateOrder(afpChain, ca1);
					
					long endTime = System.nanoTime();
					
					durationSingle = (endTime - startTime);
				}
				else if (j==2){
					//Perform the alignment
					CESymmParameters params = new CESymmParameters();
					params.setRefineMethod(RefineMethod.NOT_REFINED);
					
					afpChain = ceSymm.align(ca1, ca2, params);
					
					long endTime = System.nanoTime();
					
					durationNoRefine = (endTime - startTime);
				}
				else {
					//Perform the alignment
					CESymmParameters params = new CESymmParameters();
					params.setRefineMethod(RefineMethod.MULTIPLE);
					params.setMaxNrIterationsForOptimization(3);
					
					afpChain = ceSymm.align(ca1, ca2);
					
					orderM = afpChain.getBlockNum();
					
					long endTime = System.nanoTime();
					
					durationMultipleOpt = (endTime - startTime);
					length = ca1.length;
				}
			}
			writer.append(name+","+orders[i]+","+orderM+","+orderS+","+length+","+durationMultiple+","+durationSingle+","+durationNoRefine+","+durationMultipleOpt+"\n");
		}
		
	    writer.flush();
	    writer.close();
	}
	
	public static void main(String[] args) throws StructureException, IOException, OrderDetectionFailedException{
		
		String[] names = {"2F9H.A", "1SQU.A", "3HDP", "2AFG.A", "4DOU", "1HCE", "1TIE", "4I4Q", "1GEN", "1HXN", "1G61.A","1TL2.A","2JAJ.A", "1U6D", "1JOF.A", "1JTD.B", "1A12.A", "2I5I.A", "1GOT.B", "1TIM.A", "1VZW", "1NSJ"};
		int[] realOrders = {2,2,2,3,3,3,3,3,4,4,5,5,5,6,7,7,7,7,7,8,8,8};
		analyzeRunningTime(names, realOrders, "/home/scratch/stats/complexity.csv");
	}
}