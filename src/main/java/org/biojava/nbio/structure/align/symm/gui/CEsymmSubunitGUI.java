package org.biojava.nbio.structure.align.symm.gui;


import java.awt.Color;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.refine.MultipleAlignRefiner;
import org.biojava.nbio.structure.align.symm.refine.SingleAlignRefinement;
import org.biojava.nbio.structure.align.symm.subunit.SubunitTools;
import org.jcolorbrewer.ColorBrewer;

/**
 * Displays the alignment of the protein with the refined optimal alignment with the subunit blocks.
 * 
 * Tried and worked for:
 * 		- order 2: 2J5A, 3HDP, 4HHB, 1SQU.A, 2F9H.A, 3DDV.A, 4FI3.F, 1H9M.A, 1MP9.A, 1AIN, 1VYM.A, 4HPL.A, 1UBI, 1GUA.B, 1BA2.A
 * 		- order 3: 4DOU, 1VYM, 2AFG.A, 1HCE, 1TIE, 4I4Q
 *   	- order 4: 1GEN, 1HXN
 *  	- order 5: 1G61.A, 1TL2.A
 *  	- order 6: 1U6D, 
 *      - order 7: 1JOF.A, 1JTD.B, 1K3I.A, 2I5I.A, 1JV2.A, 1GOT.B, 1A12.A
 *      - order 8: 1TIM.A, 1VZW, 1NSJ
 *      - helical: 1B3U.A, 1EZG.A, 1DFJ.I, 1AWC.B, 1D0B.A
 *      - unknown: 1WD3, 1Z7X, 1DCE
 *  	
 * Did not work for:  1VYM.A (buggy rotation axis)
 *                    
 *                  
 * BUGS:   1*- For the 1JTD.B structure, the blackout is not done properly in the last alignment and the alignment is made
 *            over a black area, that is why the order of symmetry is incorrectly determined to 8 instead of 7. To reproduce
 *            the error use CEsymmColorGUI, which displays the dotplot of the last alignment (same with 1TIM.A). Solved with the
 *            matrix listener in the OrigM align method in the CeSymm class, that maintains the black regions during optimization.
 *         2*- The 3D alignment deletes some regions in the rotated (second) protein, which may be the unaligned regions between
 *            the subunits. It might be a problem with AlignmentTools.updateSuperposition. Some examples are: 1G61.A, 3HDP. Solved
 *            by clonning the AFPChain instead of initializing it from 0, the method needs the original one.
 *         3*- The information of gaps and RMSD in the Sequence Alignment Display is incorrect (it is set to 0). updateSuperposition 
 *            has to be changed (the dummy code), in order to update the values.
 *         4*- In the molecule 1G61.A the helices are not aligned although they seem to be symmetric, consider a less restrictive
 *            last step to select them also as symmetric. Solved with a last step to select the groups that do not form cycles.
 *         5*- The alignment panel does not color correctly the alignment, because it considers a " " in the symb alignment not as 
 *            a mismatch (when the alignment is not FatCat only) but rather as an aligned pair and it colors them (modification in 
 *            the biojava code DisplayAFP line 76), although they are not really aligned. Possible implications for biojava?
 *         6*- When the subunits are colored in the 3D structure, in some structures the color does not change between two subunits,
 *            they are either all blue or all green. This happens with the 1JTD.B structure (maybe red is not added properly).
 *         7- Rotation axis is not displayed with an issue with the jmol window when evaluating the string. Examples: 3DDV.A
 *         8*- The subunit selection seems to be very restrictive for proteins of higher order of symmetry. One solution could be
 *            to consider, if there are not <order> cycles, cycles of smaller size and establish (or group) the subunits by pairwise
 *            (or more) similarity groups. Different approach for order 6-8. Examples: 1VZW, 1TIM.A
 *         9- For small proteins, the TM score is very low even for the first alignment, which results to incorrectly determine the 
 *            order (higher than it is). This could be fixed by determining a threshold that considers the length of the protein.
 *            Examples: 1GUA.B, 1UBI.
 *        10*- From the alignment panel, an error is thrown because a position cannot be matched. Example: 1A12.A. Unknown but solved.
 *        11*- Protein 1G61.A gives some problems in the alignment: the colors are misplaced in the boundaries of the subunits and
 *            the FATCAT result is not properly shown (it shows the text alignment instead). Something to do with a double gap that
 *            might be ignore in the first subunit. Solved in getBlockNr of biojava display code.
 *        12*- The getAlign method does not consider double gaps when calculating the alignment strings and that is why some errors
 *            occur in the sequence alignment Display (color not correct and repeated residues). The problem was actually in the
 *            OptAln, because the residues were not contiguous in all the subunits. Solved by checking consistency between groups in
 *            the refinement method, there was a bug in the names of variables.
 *        13- In 1VYM structure there is a loop identified as not aligned, but it is present in the three subunits and the sequence
 *            is highly conserved, although in the 3D alignment is only seems to align well two of the three loops. In this case the
 *            subunit conditions are too restrictive.
 *         
 *         * Solved!
 *                    
 *  
 * @author Aleix Lafita
 * 
 * Last modified: 16.03.2015
 *
 */
public class CEsymmSubunitGUI {
	public static void main(String[] args) throws Exception{
		
		//String[] names = {"2F9H.A", "1SQU.A", "3HDP", "2AFG.A", "4DOU", "1GEN", "1G61.A", "1U6D", "1JOF.A", "1JTD.B", "1TL2.A", "2I5I.A", "1GOT.B", "1VZW", "1NSJ");
		String[] names = {"1VZW"};
		
		for (int i=0; i<names.length; i++){
			
			//Set the name of the protein structure to analyze
			System.out.println("Analyzing protein "+names[i]);
			AtomCache cache = new AtomCache();
			String name = names[i];

			//Parse atoms of the protein into two DS
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			System.out.println("Protein length: "+ca1.length);
			
			//Initialize a new CeSymm class and its parameters and a new alignment class
			CeSymm ceSymm = new CeSymm();
			AFPChain afpChain = new AFPChain();
			//MultipleAlignRefiner refiner = new MultipleAlignRefiner();
			//SingleAlignRefinement refiner = new SingleAlignRefinement();
			//ceSymm.setRefiner(refiner);
			
			/*//Set the maximum number of iterations for the cases where the order detection fails
			int maxNr = 8;
			CESymmParameters params = new CESymmParameters();
			params.setMaxNrAlternatives(maxNr);*/
			
			//Perform the alignment and store
			afpChain = ceSymm.align(ca1, ca2);
			
			//Set the colors (OPTIONS: Spectral (soft transition), Set1 (radical difference), Set2 (softer than 1), 
			//Set3 (soft colors), Paired (pairs of two colors)
			Color[] colors = ColorBrewer.Set1.getColorPalette(afpChain.getBlockNum());
			//Color[] colors = {Color.blue, Color.yellow, Color.cyan, Color.orange, Color.green, Color.magenta,  Color.pink, Color.red}; 
			afpChain.setBlockColors(colors);
			
			SubunitTools.displayColorSubunits(afpChain, name, ca1, ca2);
			//SubunitTools.displaySuperimposedSubunits(afpChain, name, ca1, ca2);

		}
	}
	
}