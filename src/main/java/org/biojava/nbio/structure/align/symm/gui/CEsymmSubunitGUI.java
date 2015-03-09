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
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.subunit.SubunitTools;

/**
 * Displays the alignment of the protein with the refined optimal alignment with the subunit blocks.
 * 
 * Tried and worked for: {4DOU-3, 3HDP-2, 4HHB-2, 1SQU-2, 2F9H-2, 3DDV-2, 3DDV.A-2, 4FI3.F-2, 1H9M.A-2, 1MP9.A-2, 1VYM-3}
 * Array out of bounds: 1TL2-2
 * Did not work for: {1G61.A-5 (bug in alignment, early termination), 1JTD.B-7 (buggy 8th alignment, not correct order detection), 
 *                    2FEE (takes long), 1VYM.A (buggy rotation axis)}
 *  
 * @author Aleix Lafita
 * 
 * Last modified: 09.03.2015
 *
 */
public class CEsymmSubunitGUI {
	public static void main(String[] args){

		//Set the name of the protein structure to analyze
		AtomCache cache = new AtomCache();
		String name = "1G61.A";

		try {
			
			//Parse atoms of the protein into two DS
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);
			System.out.println("Structure downloaded...");
						
			//Initialize a new CeSymm class and its parameters and a new alignment class
			CeSymm ceSymm = new CeSymm();
			AFPChain afpChain = new AFPChain();
			
			//Perform the alignment and store it in allAlignments
			afpChain = ceSymm.alignMultiple(ca1, ca2);
			afpChain.setName1(name);
			afpChain.setName2(name);
			
			//Display the alignment of the subunits, to evaluate the correctness of all alignments
			StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			
			//Set the rotation axis of the symmetry
			RotationAxis axis = new RotationAxis(afpChain);
			jmolPanel.evalString(axis.getJmolScript(ca1));
			
			/*
			jmolPanel.evalString("select *; spacefill off; wireframe off; backbone off");
			jmolPanel.evalString("select model=1.1 and (atomno >= "+ca1[intervals.get(0)].getPDBserial()+" and atomno <= "+ca1[intervals.get(1)].getPDBserial()+"); cartoon on");
			jmolPanel.evalString("select model=1.2 and (atomno >= "+ca1[intervals.get(2)].getPDBserial()+" and atomno <= "+ca1[intervals.get(3)].getPDBserial()+"); cartoon on");
			 */
			
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
}