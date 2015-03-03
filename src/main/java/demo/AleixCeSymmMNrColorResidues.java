package demo;

import java.util.ArrayList;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.subunit.SubunitTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;

/**
 * Demo for the CE-symm identifying and coloring subunits:
 * 
 * 1- With the order of symmetry given, perform multiple alignments with the blackout technique.
 * 2- Calculate the residues of each subunit by taking the residues in each interval present in all the alignments (consensus).
 * 3- Show in jmol the structure with the subunits colored differently.
 * 
 * Assumes the new align function for the CeSymm class that accepts a list of AFP alignments as input.
 * 
 * Tried and worked for: {4DOU-3, 3HDP-2, 1TL2-5, 4HHB-2, 1SQU-2, 2F9H-2, 3DDV-4, 4FI3.F-2, 1H9M.A-2, 1MP9.A-2, 1JTD.B-7, 1G61.A-5}
 * Did not work for: {2FEE (takes long), 1VYM.A (buggy rotation axis), 1VYM (takes long, buggy intervals)}
 * 
 * @author aleix
 *
 */

public class AleixCeSymmMNrColorResidues {

	public static void main(String[] args){

		//Set the name of the protein structure to analyze
		AtomCache cache = new AtomCache();
		String name = "1TL2";
		
		//Set the order of symmetry of the protein
		int order = 5;

		try {
			
			//Parse atoms of the protein into two DS
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);
			
			//List that contains all the AFP alignments
			ArrayList<AFPChain> afpAlignments= new ArrayList<AFPChain>();
			
			//Initialize a new CeSymm class and its parameters and a new alignment class
			CeSymm ceSymm = new CeSymm();
			CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
			AFPChain afpChain = new AFPChain();
			
			//Set the number of alternatives (blackouts)
			params.setMaxNrAlternatives(order);
			
			//Perform the alignment and store it in allAlignments
			afpChain = ceSymm.align(ca1, ca2, params, afpAlignments);
			afpChain.setName1(name);
			afpChain.setName2(name);
			System.out.println("There are "+afpAlignments.size()+" alignments.");
			for (AFPChain afp:afpAlignments){
				System.out.println("Alignment has length of "+afp.getOptLength());
			}
			
			///Use the method defined below to extract the subunit residues from the alignments
			ArrayList<ArrayList<Integer>> subunits = SubunitTools.extractSubunits(ca1, afpAlignments);
			System.out.println("Number of subunits: "+subunits.size());
			
			//Display the protein structure in jmol
			Structure structure = StructureIO.getStructure(name);
			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(structure);
			
			//Set the starting format of the protein
			jmol.evalString("select *; color lightgrey;");
			jmol.evalString("select *; spacefill off; wireframe off; cartoon on; ");
			jmol.evalString("select ligands; cartoon off;");
			
			//Loop through every subunit identified and through all its residues and color them differently
			for (int k=0; k<order; k++){
				
				//Get the total number of residues of the kth subunit
				int n = subunits.get(k).size();
				System.out.println("Number of residues of subunit "+(k+1)+" is "+n);
				String[] colors = {"cornflowerblue","green","mediumpurple","gold","indianred","lightskyblue","lightsalmon","lightcoral"};
				int residue = 0;
				String chain = "";
				
				for (int j=0; j<n; j++){
					residue = ca1[subunits.get(k).get(j)].getGroup().getResidueNumber().getSeqNum();
					chain = ca1[subunits.get(k).get(j)].getGroup().getResidueNumber().getChainId();
					jmol.evalString("select "+residue+":"+chain+"; color "+colors[k]);
				}
			}
						
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(afpChain);
			jmol.evalString(axis.getJmolScript(ca1));
			
			//Also display the last alignment of the subunits, to evaluate the correctness of all alignments
			StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			jmolPanel.evalString(axis.getJmolScript(ca1));
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * Method that converts a set of alignments (corresponding to the symmetry rotations) into the different subunits of the protein.
	 * 
	 * INPUT: a list of AFPChain.
	 * OUTPUT: a list of residue blocks representing each subunit.
	 * ALGORITHM: It extracts the subunit intervals from the alignments and then includes as part of every subunit only the 
	 *            residues that are present in all the alignments that contain that subunit.
	 */
	public static ArrayList<ArrayList<Integer>> extractSubunits(Atom[] ca1, ArrayList<AFPChain> allAlignments){
	
	
		//Convert the AFP alignments into a list of lists of alignment blocs called afpResidues
		ArrayList<ArrayList<Integer>> afpResidues = new ArrayList<ArrayList<Integer>>();
		int order = allAlignments.size()+1;
		
		//Get all the residue numbers for each alignment and store them as lists of Integers inside afpResidues
		//Repeat two times (y) because the afpChain divides the protein into two alignment groups (consider the 4?)
		for (int y=0; y<2; y++){
			for (int k=0; k<(order-1); k++){
				
				//Get the total number of residues of the kth subunit from the kth alignment
				int n = allAlignments.get(k).getOptAln()[0][y].length;
				ArrayList<Integer> residues = new ArrayList<Integer>();
				
				//Append every residue to a list of residues and add it to afpResidues
				for (int j=0; j<n; j++){
					residues.add(allAlignments.get(k).getOptAln()[0][y][j]);
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
					if (k==0)
						afpResidues.add(residues2);
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
				if (k==0)
					afpResidues.add(residues);
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
		System.out.println("Number of afpResidues: "+afpResidues.size());
		
		//Overlapping intervals problem: we have a set of 2*n intervals that overlap and define n subunits
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
}
