package org.biojava.nbio.structure.align.symm.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.symm.axis.SymmetryAxes;

/**
 * Action Listener for the symmetry menu.
 * Trigger various symmetry analysis.
 * 
 * @author Aleix Lafita
 *
 */
public class SymmetryListener implements ActionListener{

	private MultipleAlignmentJmol jmol;
	private MultipleAlignment msa;
	private SymmetryAxes axes;

	public SymmetryListener(MultipleAlignmentJmol jmol, SymmetryAxes axes) {
		this.jmol = jmol;
		this.msa = jmol.getMultipleAlignment();
		this.axes = axes;
	}

	@Override
	public void actionPerformed(ActionEvent ae) {

		String cmd = ae.getActionCommand();
		if (cmd.equals("Subunit Superposition")){
			if (msa == null) {
				System.err.println("Currently not displaying a symmetry!");
				return;
			}
			try {
				SymmetryDisplay.displaySubunits(msa);
			} catch (StructureException e) {
				e.printStackTrace();
			}

		} else if (cmd.equals("Multiple Structure Alignment")){
			if (msa == null) {
				System.err.println("Currently not displaying a symmetry!");
				return;
			}
			try {
				SymmetryDisplay.displayFull(msa);
			} catch (StructureException e) {
				e.printStackTrace();
			}

		} else if (cmd.equals("Point Group Symmetry")){
			if (msa != null) {
				String script = SymmetryDisplay.printPointGroupAxes(msa);
				jmol.evalString(script);
				return;
			}
			
		} else if (cmd.equals("Show Symmetry Axes")){
			if (axes != null) {
				String script = SymmetryDisplay.printSymmetryAxes(msa, axes);
				jmol.evalString(script);
				return;
			} else System.err.println("No axes for this symmetry");
			
		} else if (cmd.equals("New Symmetry Analysis")){
			SymmetryGui.getInstance();
		}
	}
	
}