/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public License.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2015-02-27
 *
 */
package org.biojava.nbio.structure.align.symm.gui;


import java.util.ArrayList;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.subunit.SubunitTools;

/**
 * Displays a two subunit alignment of the protein by simply displaying the alignment and only coloring one aligned part.
 * TODO: identify the subunits that overlap in the displayed alignment to color them only.
 *  
 * @author Aleix Lafita
 *
 */
public class CEsymmSubunitGUI {
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
			
			//Set the number of alternatives (blackouts)
			params.setMaxNrAlternatives(order);
			
			//Perform the alignment and store it in allAlignments
			afpChain = ceSymm.align(ca1, ca2, params, afpAlignments);
			afpChain.setName1(name);
			afpChain.setName2(name);
			
			//Use the method defined below to extract the subunit residues from the alignments
			ArrayList<ArrayList<Integer>> afpResidues = SubunitTools.processMultipleAFP(afpAlignments);
			ArrayList<Integer> intervals = SubunitTools.calculateIntervals(ca1, afpResidues, order);
			
			//Also display the last alignment of the subunits, to evaluate the correctness of all alignments
			StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afpAlignments.get(0), ca1, ca2);
			
			jmolPanel.evalString("select *; spacefill off; wireframe off; backbone off");
			jmolPanel.evalString("select model=1.1 and (atomno >= "+ca1[intervals.get(0)].getPDBserial()+" and atomno <= "+ca1[intervals.get(1)].getPDBserial()+"); cartoon on");
			jmolPanel.evalString("select model=1.2 and (atomno >= "+ca1[intervals.get(2)].getPDBserial()+" and atomno <= "+ca1[intervals.get(3)].getPDBserial()+"); cartoon on");
		
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
}