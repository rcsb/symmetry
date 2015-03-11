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
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;

/**
 * Displays the alignment of the protein with the refined optimal alignment with the subunit blocks.
 * 
 * Tried and worked for:
 * 		- order 2: 3HDP, 4HHB, 1SQU, 2F9H, 3DDV, 3DDV.A, 4FI3.F, 1H9M.A, 1MP9.A, 1TL2, 1AIN
 * 		- order 3: 4DOU, 1VYM, 2AFG.A, 1CUN
 *   	- order 4: 
 *  	- order 5: 1G61.A
 *  	- order 6: 1U6D
 *      - order 7: 1JOF.A, 1JTD.B, 
 *      - helical: 1B3U.A, 1EZG.A, 1DFJ.I, 1AWC.B
 *      - unknown: 1WD3, 1Z7X, 1DCE,
 *  	
 * Did not work for:  1TIM.A (symmetric?)
 * 					  2FEE (takes long)
 *                    1VYM.A (buggy rotation axis)
 *                    1DCE (takes long)
 *                    1QIU (takes long)
 *                    1TIM.A (incorrect order detection, order 7)
 *                    1VZW (incorrect order detection, 4 instead of 8)
 *                  
 * BUGS:   1- For the 1JTD.B structure, the blackout is not done properly in the last alignment and the alignment is made
 *            over a black area, that is why the order of symmetry is incorrectly determined to 8 instead of 7. To reproduce
 *            the error use CEsymmColorGUI, which displays the dotplot of the last alignment (same with 1TIM.A).
 *         2*- The 3D alignment deletes some regions in the rotated (second) protein, which may be the unaligned regions between
 *            the subunits. It might be a problem with AlignmentTools.updateSuperposition. Some examples are: 1G61.A, 3HDP. Solved
 *            by clonning the AFPChain instead of initializing it from 0, the method needs the original one.
 *         3- The information of gaps and RMSD in the Sequence Alignment Display is incorrect (it is set to 0). updateSuperposition 
 *            has to be changed (the dummy code), in order to update the values.
 *         4*- In the molecule 1G61.A the helices are not aligned although they seem to be symmetric, consider a less restrictive
 *            last step to select them also as symmetric. Solved with a last step to select the groups that do not form cycles.
 *         5- The alignment panel does not color correctly the alignment, because it considers a " " in the symb alignment not as 
 *            a mismatch (when the alignment is not FatCat only) but rather as an aligned pair and it colors them (modification in 
 *            the biojava code DisplayAFP line 76), although they are not really aligned. Possible implications for biojava? 
 *         
 *         * Solved!
 *                    
 *  
 * @author Aleix Lafita
 * 
 * Last modified: 10.03.2015
 *
 */
public class CEsymmSubunitGUI {
	public static void main(String[] args) throws Exception{
		
		//String[] names = ["4DOU", "1G61.A", "1U6D", "1JOF.A"];

		//Set the name of the protein structure to analyze
		AtomCache cache = new AtomCache();
		String name = "4DOU";

		try {
			
			//Parse atoms of the protein into two DS
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);
						
			//Initialize a new CeSymm class and its parameters and a new alignment class
			CeSymm ceSymm = new CeSymm();
			AFPChain afpChain = new AFPChain();
			
			//Set the maximum number of iterations for the cases where the order detection fails
			int order = 8;
			CESymmParameters params = new CESymmParameters();
			params.setMaxNrAlternatives(order);
			
			//Perform the alignment and store
			afpChain = ceSymm.alignMultiple(ca1, ca2, params);
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