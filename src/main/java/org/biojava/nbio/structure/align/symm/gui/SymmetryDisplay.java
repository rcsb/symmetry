package org.biojava.nbio.structure.align.symm.gui;

import java.awt.event.KeyEvent;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.symm.axis.SymmetryAxes;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.symmetry.analysis.CalcBioAssemblySymmetry;
import org.biojava.nbio.structure.symmetry.core.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorPointGroup;
import org.biojava.nbio.structure.utils.SymmetryTools;

/**
 * Class that provides visualizations methods for symmetry
 * alignments. Call the display() method for the default 
 * visualization of symmetry.
 * 
 * @author Aleix Lafita
 * 
 */
public class SymmetryDisplay {

	/**
	 * Displays a multiple alignment of the symmetry subunits.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol displaySubunits(MultipleAlignment msa) 
			throws StructureException {

		MultipleAlignment subunits = SymmetryTools.toSubunitAlignment(msa);
		return MultipleAlignmentDisplay.display(subunits);
	}

	/**
	 * Displays a multiple alignment of the whole structure transformations
	 * colored by blocks, corresponding to the subunits.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol displayFull(MultipleAlignment msa) 
			throws StructureException {

		MultipleAlignment full = SymmetryTools.toFullAlignment(msa);

		MultipleAlignmentJmol jmol = MultipleAlignmentDisplay.display(full);
		jmol.setColorByBlocks(true);
		
		return jmol;
	}
	
	/**
	 * Displays a single structure in a cartoon representation with each
	 * symmetric subunit colored differently.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @param axes symmetry axes
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol display(MultipleAlignment msa,
			SymmetryAxes axes) throws StructureException {
		
		MultipleAlignmentJmol jmol = SymmetryDisplay.displayFull(msa);
		
		//Send some commands for a nicer view
		jmol.evalString("select *; backbone off; cartoon on; model 1;");
		addSymmetryMenu(jmol, axes);
		
		//Show all the axes in the initial view
		if (axes!=null) jmol.evalString(printSymmetryAxes(msa, axes));
		jmol.evalString(printPointGroupAxes(msa));
		
		return jmol;
	}
	
	/**
	 * Displays a single structure in a cartoon representation with each
	 * symmetric subunit colored differently.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol display(MultipleAlignment msa)
			throws StructureException {
		return display(msa, null);
	}
	
	/**
	 * Adds a Symmetry menu to the Jmol display, so that further symmetry
	 * analysis can be triggered.
	 * 
	 * @param jmol parent jmol
	 * @param axes symmetry axes
	 */
	private static void addSymmetryMenu(MultipleAlignmentJmol jmol, 
			SymmetryAxes axes){
		
		JMenuBar menubar = jmol.getFrame().getJMenuBar();
		
		JMenu symm = new JMenu("Symmetry");
		symm.setMnemonic(KeyEvent.VK_S);
		
		SymmetryListener li = new SymmetryListener(jmol, axes);
		
		JMenuItem subunits = new JMenuItem("Subunit Superposition");
		subunits.addActionListener(li);
		symm.add(subunits);
		
		JMenuItem multiple = new JMenuItem("Multiple Structure Alignment");
		multiple.addActionListener(li);
		symm.add(multiple);
		
		JMenuItem pg = new JMenuItem("Point Group Symmetry");
		pg.addActionListener(li);
		symm.add(pg);
		
		JMenuItem ax = new JMenuItem("Show Symmetry Axes");
		ax.addActionListener(li);
		symm.add(ax);
		
		JMenuItem news = new JMenuItem("New Symmetry Analysis");
		news.addActionListener(li);
		symm.add(news);

		menubar.add(symm, 3);
		jmol.getFrame().pack();
	}
	
	public static String printSymmetryAxes(MultipleAlignment msa, 
			SymmetryAxes axes) {

		int id = 0;
		String script = "draw axes* off; draw poly* off;";
		Atom[] atoms = msa.getEnsemble().getAtomArrays().get(0);
		
		for (Matrix4d axis : axes.getAxes()) {
			RotationAxis rot = new RotationAxis(axis);
			script += rot.getJmolScript(atoms, id);
			id++;
		}
		return script;
	}

	public static String printPointGroupAxes(MultipleAlignment symm){

		//Split the symmetric units into different chains
		Structure subunits = SymmetryTools.toQuaternary(symm);

		//Quaternary Symmetry Detection
		QuatSymmetryParameters param = new QuatSymmetryParameters();

		CalcBioAssemblySymmetry calc = 
				new CalcBioAssemblySymmetry(subunits, param);

		QuatSymmetryDetector detector = calc.orient();
		List<QuatSymmetryResults> globalResults = detector.getGlobalSymmetry();

		AxisAligner axes = AxisAligner.getInstance(globalResults.get(0));

		//Draw the axis as in the quaternary symmetry
		JmolSymmetryScriptGenerator scriptGenerator = 
				JmolSymmetryScriptGeneratorPointGroup.getInstance(axes, "g");

		String script = "set defaultStructureDSSP true; "
				+ "set measurementUnits ANGSTROMS;  select all;  "
				+ "spacefill off; wireframe off;"
				+ "set antialiasDisplay true; autobond=false; ";

		script += scriptGenerator.getOrientationWithZoom(0);
		script += scriptGenerator.drawPolyhedron();
		script += scriptGenerator.drawAxes();
		script += "draw axes* on; draw poly* on; ";

		return script;
	}

}
