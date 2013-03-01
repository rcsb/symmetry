/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
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
 * Created on 2010-01-21
 *
 */
package org.biojava3.structure.align.symm.gui;


import javax.swing.JOptionPane;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava3.structure.align.symm.CeSymm;

/**
 * Prompts the user for a structure, then displays the CE-Symm symmetry.
 * 
 * Further alignments can be made through the menus, but they will not include
 * rotation axis graphics.
 * 
 * @author Spencer Bliven
 *
 */
public class CEsymmGUI {
	private static final long serialVersionUID = -1124973530190871875L;

	public static void main(String[] args) {
		//Add CeSymm to the top of the algorithm list
		StructureAlignment[] algorithms = StructureAlignmentFactory.getAllAlgorithms();
		StructureAlignmentFactory.clearAlgorithms();
		StructureAlignmentFactory.addAlgorithm(new CeSymm());
		for(StructureAlignment alg: algorithms) {
			StructureAlignmentFactory.addAlgorithm(alg);
		}

		//Get the pdb name from the user
		String pdb = null;
		AtomCache cache = new AtomCache();
		Atom[] ca1 = null;
		Atom[] ca2 = null;

		final String default_prompt = "Input PDB ID or domain name.\n\nExamples: 3C1G, 3JUT.A, d1bwua_";
		String prompt=default_prompt;
		while(ca1 == null) {
			pdb = (String)JOptionPane.showInputDialog(
					null,
					prompt,
					"CE-Symm",
					JOptionPane.PLAIN_MESSAGE);
//			pdb = "1VR8";
//			pdb = "3C1G";
//			pdb = "3JUT.A";
//			pdb = "d1bwua_";
			
			if(pdb == null) return; //User cancel
			else if(pdb.length()==0) continue; // Empty
			
			try {
				ca1 = cache.getAtoms(pdb);
				ca2 = cache.getAtoms(pdb);
			} catch (Exception e) {
				String error = e.getMessage();
				if(error == null) {
					error = "Error";
				}
				prompt = String.format("%s%n%s", error,default_prompt);
				pdb = null;
			}

		}
		
		// Perform the CESymm alignment
		try {
			StructureAlignment cesymm = StructureAlignmentFactory.getAlgorithm(CeSymm.algorithmName);

			AFPChain afp = cesymm.align(ca1, ca2);
			afp.setName1(pdb);
			afp.setName2(pdb);
			
			
			RotationAxis axis = new RotationAxis(afp);
			StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afp, ca1, ca2);

			String cmd = axis.getJmolScript(ca1);
			jmolPanel.evalString(cmd);
		} catch(StructureException e) {
			e.printStackTrace();
		}


	}
}