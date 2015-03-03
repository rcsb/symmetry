package org.biojava.nbio.structure.align.symm.gui;

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
 * CE-symm identifying and coloring subunits:
 * 
 * 1- With the order of symmetry given, perform multiple alignments with the blackout technique.
 * 2- Calculate the residues of each subunit by taking the residues in each interval present in all the alignments (consensus).
 * 3- Show in jmol the structure with the subunits colored differently.
 * 
 * Assumes the new align function for the CeSymm class that accepts a list of AFP alignments as input.
 * 
 * Tried and worked for: {4DOU-3, 3HDP-2, 1TL2-5, 4HHB-2, 1SQU-2, 2F9H-2, 3DDV-4, 4FI3.F-2, 1H9M.A-2, 1MP9.A-2, 1JTD.B-7, 1G61.A-5}
 * Did not work for: {2FEE (takes long), 1VYM.A (buggy rotation axis), 1VYM-6 (buggy intervals)}
 * 
 * @author Aleix Lafita
 * 
 * Last modified: 03.03.2015
 */

public class CEsymmColorGUI {

	public static void main(String[] args){

		//Set the name of the protein structure to analyze
		AtomCache cache = new AtomCache();
		String name = "4dou";
		
		//Set the order of symmetry of the protein
		int order = 3;

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
			
			//Set the number of alternatives (blackouts) for this iteration
			params.setMaxNrAlternatives(order);
			
			//Perform the alignment and store it in allAlignments
			afpChain = ceSymm.align(ca1, ca2, params, afpAlignments);
			afpChain.setName1(name);
			afpChain.setName2(name);
			System.out.println("There are "+afpAlignments.size()+" alignments.");
			for (AFPChain afp:afpAlignments){
				System.out.println("Alignment has length of "+afp.getOptLength());
			}
			
			//Use the method defined above to extract the subunit residues from the alignments
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
}
