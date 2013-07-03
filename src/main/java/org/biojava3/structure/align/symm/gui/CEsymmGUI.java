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
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.scop.ScopFactory;
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
		
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75B);
		
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
		Atom[] caInter = null;
		
		final String default_prompt = "Input PDB ID or domain name.\n\nExamples: 3C1G, 3JUT.A, d1bwua_";
		String prompt=default_prompt;
		while(ca1 == null) {
//			pdb = "1VR8";
//			pdb = "3C1G";
//			pdb = "3JUT.A";
//			pdb = "d1bwua_";
//			pdb = "d2biba1"; // screw symmetry
//			pdb = "d1lghb_"; // negative screw symmetry
//			pdb = "d1v3wa_"; // purely translational
//			pdb = "d1yqha1"; // buggy case
//			pdb = "d1fwka2";
//			pdb = "d1ewfa1";

			pdb = (String)JOptionPane.showInputDialog(
					null,
					prompt,
					"CE-Symm",
					JOptionPane.PLAIN_MESSAGE);
			
			if(pdb == null) return; //User cancel
			else if(pdb.length()==0) continue; // Empty
			
			try {
				ca1 = cache.getAtoms(pdb);
				ca2 = cache.getAtoms(pdb);
				caInter = cache.getAtoms(pdb);
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
			
			Matrix mat = afp.getBlockRotationMatrix()[0];
			Atom shift = afp.getBlockShiftVector()[0];
			for( Atom atom:caInter) {
				Calc.rotate(atom, mat);
			}
			
			RotationAxis axis = new RotationAxis(afp);
			StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afp, ca1, ca2);
			
			String cmd = axis.getJmolScript(ca1);
			jmolPanel.evalString(cmd);
			
			
			System.out.println("Theta="+axis.getAngle());
			/*
			jmolPanel.evalString("set perspectiveDepth 0;"); //orthoscopic 
			StringBuilder debugCmds = new StringBuilder();
			debugCmds.append(String.format(
					"draw ID ca1 arrow {%f,%f,%f} {%f,%f,%f} width 0.5 \">ca1\";%n",
					ca1[0].getX(), ca1[0].getY(),ca1[0].getZ(),
					ca1[ca1.length-1].getX(), ca1[ca1.length-1].getY(),ca1[ca1.length-1].getZ()));
			debugCmds.append(String.format(
					"draw ID ca1b arrow {%f,%f,%f} {%f,%f,%f} width 0.5 \">ca1b\";%n",
					ca1[0].getX(), ca1[0].getY(),ca1[0].getZ(),
					ca1[16].getX(), ca1[16].getY(),ca1[16].getZ()));
			debugCmds.append(String.format(
					"draw ID ca2 arrow {%f,%f,%f} {%f,%f,%f} width 0.5 \">ca2\";%n",
					ca2[0].getX(), ca2[0].getY(),ca2[0].getZ(),
					ca2[ca2.length-1].getX(), ca2[ca2.length-1].getY(),ca2[ca2.length-1].getZ()));
			debugCmds.append(String.format(
					"draw ID ca2b arrow {%f,%f,%f} {%f,%f,%f} width 0.5 \">ca2b\";%n",
					ca2[0].getX(), ca2[0].getY(),ca2[0].getZ(),
					ca2[16].getX(), ca2[16].getY(),ca2[16].getZ()));
			debugCmds.append(String.format(
					"draw ID caI arrow {%f,%f,%f} {%f,%f,%f} width 0.5 \">caI\";%n",
					caInter[0].getX(), caInter[0].getY(),caInter[0].getZ(),
					caInter[caInter.length-1].getX(), caInter[caInter.length-1].getY(),caInter[caInter.length-1].getZ()));
			debugCmds.append(String.format(
					"draw ID caIb arrow {%f,%f,%f} {%f,%f,%f} width 0.5 \">caIb\";%n",
					caInter[0].getX(), caInter[0].getY(),caInter[0].getZ(),
					caInter[16].getX(), caInter[16].getY(),caInter[16].getZ()));
			jmolPanel.evalString(debugCmds.toString());

			
			// draw intermediate vectors for debugging
			double width = .5;
			Atom s = axis.getRotationPos();
			Atom u = axis.getRotationAxis();
			jmolPanel.evalString(String.format("draw ID s VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR orange \">s\";",
					s.getX(),s.getY(),s.getZ(), width ));

			Atom perp = axis.getOtherTranslation();
			Atom screw = axis.getScrewTranslation();

			double uScale = 10;
			
			jmolPanel.evalString(String.format("draw ID u VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR orange \">u\";",
					uScale*u.getX(),uScale*u.getY(),uScale*u.getZ(), width ));

			jmolPanel.evalString(String.format("draw ID perp VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR yellow \">tPerp\";",
					perp.getX(),perp.getY(),perp.getZ(), width));
			jmolPanel.evalString(String.format("draw ID screw VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR yellow \">screw\";",
					screw.getX(),screw.getY(),screw.getZ(), width));
			
			// draw coordinate axes
			jmolPanel.evalString("draw ID x VECTOR {0,0,0} {5,0,0} WIDTH 0.5 COLOR red \">x\";");
			jmolPanel.evalString("draw ID y VECTOR {0,0,0} {0,5,0} WIDTH 0.5 COLOR green \">y\";");
			jmolPanel.evalString("draw ID z VECTOR {0,0,0} {0,0,5} WIDTH 0.5 COLOR blue \">z\";");
			*/
		} catch(StructureException e) {
			e.printStackTrace();
		}


	}
}