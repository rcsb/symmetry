package demo;

import java.util.Vector;

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

/**
 * Demo for the CE-symm coloring subunits.
 *
 */
public class AleixCeSymmMaxNr {

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
			
			//List that contains all the possible alignments 
			AFPChain[] afpAlignments= new AFPChain[order-1];
			
			//Iterate for every possible rotation number of the molecule (order) and obtain the AFPchain
			int i = 1;
			while (i<order) {
				
				CeSymm ceSymm = new CeSymm();
				//Initialize a variable for the CeSymm class parameters
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
			//Get the residue numbers (start and end) for each subunit
			int[] afpResidues = new int[order-1];
			for (int k=0; k<afpAlignments.length; k++){
				int n = afpAlignments[k].getOptAln()[0][0].length;
				afpResidues[k] = afpAlignments[k].getOptAln()[0][0][n-1];
				System.out.println(afpResidues[k]);
			}
			
			//Display the structure in jmol
			Structure structure = StructureIO.getStructure(name);
			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(structure);
			
			//Color the different chains and set style format
			jmol.evalString("select *; color chain;");
			jmol.evalString("select *; spacefill off; wireframe off; cartoon on; ");
			jmol.evalString("select ligands; cartoon off; wireframe 0.3; spacefill 0.5; color cpk;");
			int current_index = 1;
			String[] colors = {"blue", "magenta", "green", "cian"};
			for (int k=0; k<afpAlignments.length; k++){
				jmol.evalString("select "+current_index+"-"+afpResidues[k]+"; color "+colors[k]);
				current_index = afpResidues[k];
			}
						
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(afpAlignments[0]);
			jmol.evalString(axis.getJmolScript(ca1));
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
