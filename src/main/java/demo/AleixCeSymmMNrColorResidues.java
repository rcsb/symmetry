package demo;

import java.util.Vector;
import java.util.Arrays;

import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.junit.runner.manipulation.Sorter;

/**
 * Demo for the CE-symm coloring subunits.
 * With the order of symmetry given, perform multiple alignments and calculate the aligned residues of each subunit.
 * (Assumes the new align function for the CeSymm class that accepts an AFPChain array as input.)
 * Then show in jmol the structure with the subunits colored differently.
 * 
 * Tried and worked for: {4DOU, 3HDP, 4HHB, 1SQU, 2F9H, 3DDV, 1VYM (takes long), 4FI3.F}
 * Did not work for: {2FEE, 1VYM.A}
 * 
 * @author aleix
 *
 */
public class AleixCeSymmMNrColorResidues {

	@SuppressWarnings("null")
	public static void main(String[] args){

		//Set the name of the protein structure to analyze
		AtomCache cache = new AtomCache();
		String name = "4DOU";
		
		//Set the order of symmetry of the protein
		int order = 10;

		try {
			
			//Parse atoms of the protein into two DS
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);	
			
			//Array that contains all the possible alignments 
			AFPChain[] afpAlignments= new AFPChain[order-1];
			
			//Iterate for every possible rotation number of the molecule (order) and obtain the AFPchain
			int i = 1;
			while (i<order) {
				
				//Initialize a new CeSymm class and its parameters and a new alignment class
				CeSymm ceSymm = new CeSymm();
				CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
				AFPChain afpChain = new AFPChain();
				
				//Set the number of alternatives (blackouts) for this iteration
				params.setMaxNrAlternatives(i);
				
				//Perform the alignment and store it in allAlignments
				afpChain = ceSymm.align(ca1, ca2, params);
				afpChain.setName1(name);
				afpChain.setName2(name);
				System.out.println(afpChain.getOptLength());
				afpAlignments[i-1] = afpChain;
				i++;
			}
			
			//Get the residue numbers for each subunit - increment the loop to get all the residue numbers
			int[] afpResidues = new int[order]; //array containing the residue indices
			afpResidues[0] = ca1[0].getGroup().getResidueNumber().getSeqNum();
			System.out.println(ca1.length);
					
			for (int k=0; k<order-1; k++){
				int n = afpAlignments[k].getOptAln()[0][0].length;
				afpResidues[k+1] = ca1[afpAlignments[k].getOptAln()[0][0][n-1]].getGroup().getResidueNumber().getSeqNum();
				System.out.println(afpResidues[k+1]);
			}
			
			//Display the structure in jmol
			Structure structure = StructureIO.getStructure(name);
			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(structure);
			
			//Color the different chains and set style format
			jmol.evalString("select *; color [0,0,200];");
			jmol.evalString("select *; spacefill off; wireframe off; cartoon on; ");
			jmol.evalString("select ligands; cartoon off;");
			for (int k=0; k<order-1; k++){
				int colorb = (200-(k+1)*(200/(order)));
				int colorg = ((k+1)*(200/(order-1)));
				jmol.evalString("select "+afpResidues[k]+"-"+afpResidues[k+1]+"; color [0,"+colorg+","+colorb+"]");
			}
						
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(afpAlignments[0]);
			jmol.evalString(axis.getJmolScript(ca1));
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
