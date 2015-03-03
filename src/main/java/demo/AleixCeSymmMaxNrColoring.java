package demo;

import java.util.Arrays;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;

/**
 * Demo for the CE-symm coloring subunits.
 * With the order of symmetry given, perform multiple alignments and calculate the start and end atoms of each subunit.
 * Then show in jmol the structure with the subunits colored differently.
 * 
 * Tried and worked for: 4DOU, 3HDP, 4HHB, 1SQU, 2F9H, 3DDV, 1VYM (takes long), 4FI3.F
 * Did not work for: 2FEE, 1VYM.A
 * 
 * @author aleix
 *
 */
public class AleixCeSymmMaxNrColoring {

	@SuppressWarnings("null")
	public static void main(String[] args){

		//Set the name of the protein structure to analyze
		AtomCache cache = new AtomCache();
		String name = "4DOU";
		
		//Set the order of symmetry of the protein
		int order = 3;

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
			
			//Get the atom numbers (start and end) for each subunit
			int[] afpResidues = new int[order]; //array containing the atom indices
			afpResidues[0] = ca1[0].getPDBserial();
			System.out.println(ca1.length);
			System.out.println(ca2.length);
					
			for (int k=0; k<order-1; k++){
				int n = afpAlignments[k].getOptAln()[0][0].length;
				afpResidues[k+1] = ca1[afpAlignments[k].getOptAln()[0][0][n-1]].getPDBserial();
				System.out.println(afpResidues[k+1]);
			}
			//Sort the atom indices from lower to higher
			//Issue: it does not ensure correctness, because different alignments can start at the same atom
			Arrays.sort(afpResidues);
			
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
				jmol.evalString("select atomno >= "+afpResidues[k]+" and atomno <= "+afpResidues[k+1]+"; color [0,"+colorg+","+colorb+"]");
			}
						
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(afpAlignments[0]);
			jmol.evalString(axis.getJmolScript(ca1));
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
