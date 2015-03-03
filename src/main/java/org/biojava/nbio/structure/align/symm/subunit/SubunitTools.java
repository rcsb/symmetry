package org.biojava.nbio.structure.align.symm.subunit;

import java.util.ArrayList;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * Methods to process AFPChain multiple alignments and analyze the symmetrical subunits generated.
 * 
 * @author Aleix Lafita
 * 
 * Last modified: 03.03.2015
 */
public class SubunitTools {
	
	/**
	 * Method that converts a set of alignments (corresponding to the symmetry rotations) into the different subunits of the protein.
	 * 
	 * INPUT: a list of AFPChain alignments and the protein atoms.
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It extracts the subunit intervals from the alignments and then includes as part of every subunit only the 
	 *            residues that are present in all the alignments that contain that subunit range.
	 */
	public static ArrayList<ArrayList<Integer>> extractSubunits(Atom[] ca1, ArrayList<AFPChain> allAlignments){
	
		//Set the order of symmetry, process the multiple alignment and calculate the intervals
		int order = allAlignments.size()+1;
		ArrayList<ArrayList<Integer>> afpResidues = processMultipleAFP(allAlignments);
		ArrayList<Integer> intervals = calculateIntervals(ca1, afpResidues,order);
		
		//Define the residues of each subunit with the intervals (cut regions for the subunits) and only add each residue
		//if it is contained in all the alignment blocks that contain that subunit (interval)
		ArrayList<ArrayList<Integer>> subunits = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> subunit = new ArrayList<Integer>();
		//sub variable keeps track of the current subunit
		int sub = 1;
		
		for (int i=0; i<ca1.length; i++){
			
			//include determines if the residue is part of a subunit or not
			boolean include = false;
			
			for (int k=0; k<afpResidues.size(); k++){
				//If the residue is between the first and last residues of an alignment block do
				if (i>afpResidues.get(k).get(0) && i<afpResidues.get(k).get(afpResidues.get(k).size()-1)){
					//If it is included set include to true temporarily and continue with the next block
					if (afpResidues.get(k).contains(i)){
						include = true;
					}
					//If it is not included break the for loop and go to the next residue
					else {
						include = false;
						break;
					}
				}
			}
			if (include){
					subunit.add(i);
			}
			if (intervals.get(sub)==i){
				ArrayList<Integer> copySub = new ArrayList<Integer>(subunit);
				subunits.add(copySub);
				subunit.clear();
				sub++;
			}
		}
		
	return subunits;
	}

	/**
	 * Method that calculates the intervals of every subunit (cutting regions). 
	 * 
	 * INPUT: a list of alignment blocs and an array of atoms of the protein.
	 * OUTPUT: a list containing the cutting intervals (size: order+1).
	 * ALGORITHM: It takes all the starting and ending residues from the alignment blocks and clusters them by iteratively
	 *            deleting the higher of the closest pair of residues. Other approaches to cluster are being considered.
	 * 
	 */
	public static ArrayList<Integer> calculateIntervals(Atom[] ca1, ArrayList<ArrayList<Integer>> afpResidues, int order){
				
		//Overlapping intervals problem: we have a set of 4*n intervals that overlap and define n subunits
		//Try to define the n subunits from the residue sets by its starting and ending residues of each interval
		ArrayList<Integer> intervals = new ArrayList<Integer>();
		
		//Add the first and last residue numbers of the protein
		intervals.add(0);
		intervals.add(ca1.length-1);
		
		for (int k=0; k<afpResidues.size()-1; k++){
			//Add the first and last residue numbers of each alignment block
			int n = afpResidues.get(k).size();
			intervals.add(afpResidues.get(k).get(0));
			intervals.add(afpResidues.get(k).get(n-1));				
		}
		
		//Sort the interval numbers and delete the closer ones until only there are <order+1> numbers left
		intervals.sort(null);
		while (intervals.size()>order+1){
			int toBdeleted = 0;
			int minDistance = ca1.length;
			int distance = 0;
			for (int j=1; j<intervals.size(); j++){
				distance = intervals.get(j)-intervals.get(j-1);
				if (distance < minDistance){
					minDistance = distance;
					toBdeleted = j;
				}
			}
			intervals.remove(toBdeleted);
		}
		//Ensure that the last residue number is the last residue of the molecule
		intervals.remove(intervals.size()-1);
		intervals.add(ca1.length-1);
		
		//Temporal: analyze the intervals
		System.out.println("Number of interval residues: "+intervals.size());
		for (int intr:intervals){
			System.out.println(intr);
		}
		
		return intervals;
	}
	
	/**
	 * Method that calculates the intervals of every subunit (cutting regions). 
	 * 
	 * INPUT: a list of AFP alignments.
	 * OUTPUT: a list containing all the alignment blocks from the AFP alignments, ordered by their starting residue.
	 * 
	 */
	public static ArrayList<ArrayList<Integer>> processMultipleAFP(ArrayList<AFPChain> allAlignments){
	
		//Convert the AFP alignments into a list of lists of alignment blocs called afpResidues
		ArrayList<ArrayList<Integer>> afpResidues = new ArrayList<ArrayList<Integer>>();
		int order = allAlignments.size()+1;
		
		//Get all the residue numbers for each alignment and store them as lists of Integers inside afpResidues
		//Repeat four times because the afpChain divides the protein into two alignment groups for each copy
		for (int k=0; k<(order-1); k++){
			for (int y=0; y<allAlignments.get(k).getOptAln().length; y++){
				for (int z=0; z<allAlignments.get(k).getOptAln()[y].length; z++){
							
					//Get the total number of residues of the kth subunit from the kth alignment
					int n = allAlignments.get(k).getOptAln()[y][z].length;
					ArrayList<Integer> residues = new ArrayList<Integer>();
					
					//Append every residue to a list of residues and add it to afpResidues
					for (int j=0; j<n; j++){
						residues.add(allAlignments.get(k).getOptAln()[y][z][j]);
					}
					//Special case: if the last residue is smaller than the last, split the block into two sub-blocks
					if (residues.get(0)>residues.get(residues.size()-1)){
						ArrayList<Integer> residues2 = new ArrayList<Integer>();
						for (int i=1; i<residues.size(); i++){
							if (residues.get(i)<residues.get(0)){
								residues2.add(residues.get(i));
								residues.remove(i);
							}	
						}
						//Store the subunits in increasing order
						if (k==0 && y==0 && z==0){
							afpResidues.add(residues2);
						}
						else{
							//Compare the first and last residue numbers of both intervals and decide accordingly
							for (int x=0; x<k; x++){
								if(afpResidues.get(x).get(0)>residues2.get(0)){
									afpResidues.add(x, residues2);
								}
								else if(x==k-1){
									afpResidues.add(residues2);
								}
							}
						}
					}
					//Store the subunits in increasing order
					if (k==0 && y==0 && z==0){
						afpResidues.add(residues);
					}
					else{
						//Compare the first and last residue numbers of both intervals and decide accordingly
						for (int x=0; x<k; x++){
							if(afpResidues.get(x).get(0)>residues.get(0)){
								afpResidues.add(x, residues);
							}
							else if(x==k-1){
								afpResidues.add(residues);
							}
						}
					}
				}
			}
		}
		System.out.println("Number of afpResidues: "+afpResidues.size());
		return afpResidues;
	}
	
	
	/**
	 * Method that displays two superimposed subunits in jmol. In construction...
	 * 
	 * INPUT: a list of subunits and a list of AFP alignments.
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It extracts the subunit intervals from the alignments and then includes as part of every subunit only the 
	 *            residues that are present in all the alignments that contain that subunit.
	 */
	public static void displaySubunits(ArrayList<ArrayList<Integer>> subunits, ArrayList<AFPChain> allAlignments){
		
		
	}
	
	/**
	 * Method that displays the protein and colors each subunits differently in jmol. In construction...
	 * 
	 * INPUT: a list of subunits and a list of AFP alignments.
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It extracts the subunit intervals from the alignments and then includes as part of every subunit only the 
	 *            residues that are present in all the alignments that contain that subunit.
	 */
	public static void colorSubunits(ArrayList<ArrayList<Integer>> subunits, ArrayList<AFPChain> allAlignments){
		
		
	}
	
}
