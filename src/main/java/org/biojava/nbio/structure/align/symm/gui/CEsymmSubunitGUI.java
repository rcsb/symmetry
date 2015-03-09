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


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.align.symm.CeSymm;

/**
 * Displays the alignment of the protein with the refined optimal alignment with the subunit blocks.
 * 
 * Tried and worked for:
 * 		- order 2: 3HDP, 4HHB, 1SQU, 2F9H, 3DDV, 3DDV.A, 4FI3.F, 1H9M.A, 1MP9.A, 1TL2, 1AIN
 * 		- order 3: 4DOU, 1VYM, 2AFG.A, 1CUN
 *   	- order 4: 
 *  	- order 5: 1G61.A
 *  	- order 6: 
 *      - order 7: 
 *      - helical: 1B3U.A, 1EZG.A, 1DFJ.I
 *      - unknown: 1WD3, 1Z7X, 1DCE,
 *  	
 * Did not work for:  1JTD.B-7 (buggy 8th alignment, not correct order detection)
 *                    2FEE (takes long)
 *                    1VYM.A (buggy rotation axis)
 *                    1AWC
 *                    1DCE (takes long)
 *                    1QIU (takes long)
 *                    1U6D (no groups selected in the refinement, too strict)
 *                    1WD3 ((no groups selected in the refinement, too strict)
 *                    1FT2.B (refinement not sorted, getAlign fails)
 *                    1JOF.A (refinement not sorted, getAlign fails)
 *                    1VZW (incorrect order detection)
 *                    
 *  
 * @author Aleix Lafita
 * 
 * Last modified: 09.03.2015
 *
 */
public class CEsymmSubunitGUI {
	public static void main(String[] args) throws Exception{

		//Set the name of the protein structure to analyze
		AtomCache cache = new AtomCache();
		String name = "1WD3";

		try {
			
			//Parse atoms of the protein into two DS
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);
						
			//Initialize a new CeSymm class and its parameters and a new alignment class
			CeSymm ceSymm = new CeSymm();
			AFPChain afpChain = new AFPChain();
			
			//Perform the alignment and store
			afpChain = ceSymm.alignMultiple(ca1, ca2);
			afpChain.setName1(name);
			afpChain.setName2(name);
			
			//Display the AFP alignment of the subunits
			StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(afpChain);
			jmolPanel.evalString(axis.getJmolScript(ca1));		
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
}