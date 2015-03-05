package org.biojava.nbio.structure.align.symm.subunit;

import java.lang.reflect.Array;
import java.util.ArrayList;

import org.biojava.nbio.structure.align.fatcat.calc.AFPOptimizer;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.AFPTwister;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Methods to process AFPChain multiple alignments and analyze the symmetrical subunits generated.
 * 
 * @author Aleix Lafita
 * 
 * Last modified: 04.03.2015
 */
public class SubunitTools {
	
	/**
	 * Method that converts a set of alignments (corresponding to the symmetry rotations) into a list of residues 
	 * for each subunit of the protein that is consistent (equivalent starting and ending residues for subunit).
	 * 
	 * INPUT: a list of AFPChain alignments and the protein atoms.
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It extracts the subunit intervals from the alignments and then includes as part of every subunit only the 
	 *            residues that are present in all the alignment and are consistent.
	 */
	public static ArrayList<ArrayList<Integer>> extractSubunits(Atom[] ca1, ArrayList<AFPChain> allAlignments){
	
		//Process the alignments to which the subunit residues should be restricted (number higher is more restrictive)
		ArrayList<ArrayList<Integer>> afpResidues = new ArrayList<ArrayList<Integer>>();
		int number = 1; //can be in the interval [1,allAlignments.size()]
		for (int k=0; k<number; k++){
			ArrayList<Integer> residues = processAFP(allAlignments.get(k));
			afpResidues.add(residues);
		}
				
		//Calculate the intervals
		ArrayList<Integer> intervals = calculateIntervals(ca1, allAlignments);
		
		//Define the residues of each subunit with the intervals (cut regions for the subunits) and only add each residue
		//if it is contained in all the alignments
		ArrayList<ArrayList<Integer>> subunits = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> subunit = new ArrayList<Integer>();
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
				ArrayList<Integer> copySub = new ArrayList<Integer>(subunit);
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
	 * 
	 */
	public static ArrayList<Integer> calculateIntervals(Atom[] ca1, ArrayList<AFPChain> allAlignments){
				
		//Overlapping intervals problem: we have a set of 4*n intervals that overlap and define n subunits
		//Try to define the n subunits from the residue sets by its starting and ending residues of each interval
		ArrayList<Integer> intervals = new ArrayList<Integer>();
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
		ArrayList<Integer> selectedInter = new ArrayList<Integer>();
		ArrayList<Integer> alreadySeen = new ArrayList<Integer>();
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
	 * Method that processes an AFP alignment and returns a list of the residues present in both chains of the alignment.
	 * 
	 * INPUT: an AFP alignment.
	 * OUTPUT: a list containing all the aligned residues in the AFP alignment, sorted increasingly.
	 */
	public static ArrayList<Integer> processAFP(AFPChain afpChain){
		
		//Get all the residue numbers in the alignment and store them in a list
		ArrayList<Integer> residues = new ArrayList<Integer>();
		for (int y=0; y<afpChain.getOptAln().length; y++){
			
			//Get the total number of residues in the group
			int n = afpChain.getOptAln()[y][0].length;
			
			//Append every residue to a list of residues
			for (int j=0; j<n; j++){
				residues.add(afpChain.getOptAln()[y][0][j]);
			}
		}
		
		//Now do the same for the other chain, select the residues only if present in the two chains.
		ArrayList<Integer> intersectionResidues = new ArrayList<Integer>();
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
	 * INPUT: a list of subunits and the name of the protein.
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It displays the protein with jmol and colors the subunit residues using jmol commands.
	 */
	public static void colorSubunits(ArrayList<ArrayList<Integer>> subunits, String name){
		
		
	}
	
	/**
	 * It returns an AFPChain with a segmented optimal alignment, which means that it has <order of symmetry> blocks of
	 * aligned subunits.
	 * 
	 * INPUT: a list of AFP alignments and the Atom[] arrays.
	 * OUTPUT: the optimal AFPChain alignment divided into the subunits.
	 */
	public static AFPChain replaceOptAln(ArrayList<AFPChain> allAlignments, Atom[] ca1, Atom[] ca2) throws StructureException {
		
		//Extract the subunits from the set of alignments
		ArrayList<ArrayList<Integer>> subunits = extractSubunits(ca1, allAlignments);
		int order = subunits.size();
				
		//Extract the old alignment
		int[][][] oldAlgn = allAlignments.get(0).getOptAln();
		
		//Create a triple list with the optimal alignment divided into the subunits in <order> blocs
		ArrayList<ArrayList<ArrayList<Integer>>> newAlgn = new ArrayList<ArrayList<ArrayList<Integer>>>();
		
		//Loop through all the optimal alignment and subunits and select the pairs of aligned residues for each subunit
		for (int i=0; i<order; i++){
			ArrayList<Integer> chain1 = new ArrayList<Integer>();
			ArrayList<Integer> chain2 = new ArrayList<Integer>();
			for (int k=0; k<oldAlgn.length; k++){
				for (int j=0; j<oldAlgn[k][0].length; j++){
					if (subunits.get(i).get(subunits.get(i).size()-1) >= oldAlgn[k][0][j] && subunits.get(i).get(0) <= oldAlgn[k][0][j]){
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
		System.out.println("Subunits in optAlgn: "+optAlgn.length);
		for (int i=0; i<optAlgn.length; i++){
			System.out.println("Subunit length: "+optAlgn[i][0].length);
		}
		
		/*
		System.out.println("Sequence alignment block 3");
		for (int i=0; i<optAlgn[1][0].length; i++){
			System.out.println(ca1[optAlgn[1][0][i]].getGroup().getPDBName()+" - "+ca1[optAlgn[1][0][i]].getGroup().getPDBName());
		}
		*/
		
		//set everything
		AFPChain refinedAFP = (AFPChain) allAlignments.get(0).clone();
		refinedAFP.setOptLength(optLength);
		refinedAFP.setBlockSize(optLens);
		refinedAFP.setOptLen(optLens);
		refinedAFP.setOptAln(optAlgn);
		refinedAFP.setBlockNum(order);
		
		//recalculate properties: superposition, tm-score, etc
		Atom[] ca2clone = StructureTools.cloneCAArray(ca2); // don't modify ca2 positions
		AlignmentTools.updateSuperposition(refinedAFP, ca1, ca2clone);
		
		//Try to fix the sequence alignment panel bug (it does not correspond to the 3D representation - change block)
		refinedAFP.setBlockSize(optLens);
		
		//It re-does the sequence alignment strings from the OptAlgn information only (why does not work..?)
		refinedAFP.setAlnsymb(null);
		AFPAlignmentDisplay.getAlign(refinedAFP, ca1, ca2, true);
		
		return refinedAFP;
	}
}