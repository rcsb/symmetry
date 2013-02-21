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


import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.gui.RotationAxis;
import org.biojava3.structure.align.symm.CeSymm;

/**
 * 
 * @author Andreas Prlic
 *
 */
public class CEsymmGUI extends JFrame{

//	public CEsymmGUI() {
//		super("CE-Symm");
//		
//		JPanel main = createMainPanel();
//		getContentPane().add(main);
//		setDefaultCloseOperation(EXIT_ON_CLOSE);
//		pack();
//		setVisible(true);
//	}
//	private JPanel createMainPanel() {
//		JPanel main = new JPanel();
//		main.setLayout(new BorderLayout());
//		JPanel sub = new JPanel();
//		sub.add(new JLabel("PDB ID:"));
//		
//		
//		
//		main.add(sub,BorderLayout.CENTER);
//		
//		JButton submit = new JButton(new AlignAction("Get Symmetry"));
//		main.add(submit,BorderLayout.SOUTH);
//		
//		return main;
//	}
	
	private static class AlignAction extends AbstractAction {
		public AlignAction(String name) {
			super(name);
		}
		@Override
		public void actionPerformed(ActionEvent e) {

			System.out.println("Aligning ");
		}
	}
	
	public static void main(String[] args) {
	
		
		//Add CeSymm to the top of the algorithm list
		StructureAlignment[] algorithms = StructureAlignmentFactory.getAllAlgorithms();
		StructureAlignmentFactory.clearAlgorithms();
		StructureAlignmentFactory.addAlgorithm(new CeSymm());
		for(StructureAlignment alg: algorithms) {
			StructureAlignmentFactory.addAlgorithm(alg);
		}

		String pdb = (String)JOptionPane.showInputDialog(
		                    null,
		                    "PDB ID:",
		                    "CE-Symm",
		                    JOptionPane.PLAIN_MESSAGE);

		//If a string was returned, say so.
		if ((pdb != null) && (pdb.length() > 0)) {
		    try {
		    	AtomCache cache = new AtomCache();
				Atom[] ca1 = cache.getAtoms(pdb);
				Atom[] ca2 = cache.getAtoms(pdb);
				
				CeSymm ce = new CeSymm();
				
				AFPChain afp = ce.align(ca1, ca2);
				
				RotationAxis axis = new RotationAxis(afp);
				StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afp, ca1, ca2);
				
				axis.displayRotationAxis(jmolPanel, ca1);
				
				
			} catch (IOException e) {
				e.printStackTrace();
			} catch (StructureException e) {
				e.printStackTrace();
			}
		}

		
//		new CEsymmGUI();
	}
}