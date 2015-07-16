package org.biojava.nbio.structure.align.symm.axis;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymmIterative;
import org.biojava.nbio.structure.align.symm.CeSymmRecursive;
import org.biojava.nbio.structure.align.symm.ChainSorter;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.symmetry.analysis.CalcBioAssemblySymmetry;
import org.biojava.nbio.structure.symmetry.core.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorPointGroup;
import org.biojava.nbio.structure.utils.SymmetryTools;

/**
 * This class determines all symmetry axis present in a symmetry
 * alignment of the subunits. 
 * It uses the code to determine quaternary symmetry in biojava
 * as its basis.<p>
 * The algorithm splits the symmetric units of a structure into 
 * different chains first, in order to share the code for axis 
 * determination with the one for the quaternary symmetry detection.
 * 
 * @author Aleix Lafita
 *
 */
public class InternalSymmetryAxes {

	/**
	 * Calculates all symmetry axis from a symmetry alignment
	 * resulting from a CeSymm analysis.
	 * 
	 * @param symm MultipleAlignment
	 * @return axis aligner of the symmetry
	 */
	public static AxisAligner calculateAxes(MultipleAlignment symm) {

		//Split the symmetric units into different chains
		Structure subunits = SymmetryTools.toQuaternary(symm);

		//Quaternary Symmetry Detection
		QuatSymmetryParameters param = new QuatSymmetryParameters();
		param.setOnTheFly(true);
		param.setVerbose(true);
		param.setSequencePseudoSymmetryThreshold(0.0);
		param.setMinimumSequenceLengthFraction(0.0);
		param.setAlignmentFractionThreshold(0.0);
		param.setAbsoluteMinimumSequenceLength(10);
		param.setAngleThreshold(10.0);
		param.setRmsdThreshold(20.0);
		param.setMinimumSequenceLength(10);

		CalcBioAssemblySymmetry calc = 
				new CalcBioAssemblySymmetry(subunits, param);

		QuatSymmetryDetector detector = calc.orient();
		List<QuatSymmetryResults> globalResults = detector.getGlobalSymmetry();

		AxisAligner aligner = AxisAligner.getInstance(globalResults.get(0));

		return aligner;
	}

	public static void main(String[] args) throws Exception {

		//More than one symmetry axis: 4gcr, 1ppr.O, 1vym.A, 1yox.A
		//Domain swapping: 1g6s
		//Internal+quaternary: 1VYM, 1f9z, 1YOX_A:,B:,C:, 1mmi
		//Structures that have different symmetry thresholds: 1vzw
		//Dihedral structures: 4hhb, 1iy9, 2ehz,
		String name = "4gcr";

		AtomCache cache = new AtomCache();
		Atom[] atoms = ChainSorter.cyclicSorter(cache.getStructure(name));

		CESymmParameters params = new CESymmParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setOptimization(true);

		//CeSymmRecursive recurser = new CeSymmRecursive(params);
		CeSymmIterative recurser = new CeSymmIterative(params);
		MultipleAlignment msa = recurser.execute(atoms);

		AxisAligner axis = InternalSymmetryAxes.calculateAxes(msa);

		//Draw the axis as in the quaternary symmetry
		JmolSymmetryScriptGenerator scriptGenerator = 
				JmolSymmetryScriptGeneratorPointGroup.getInstance(axis, "g");

		SymmetryJmol jmol = new SymmetryJmol(msa, null);
		
		String script = "set defaultStructureDSSP true; "
				+ "set measurementUnits ANGSTROMS;  select all;  "
				+ "spacefill off; wireframe off; model 0;"
				+ "set antialiasDisplay true; autobond=false; ";
		
		script += scriptGenerator.getOrientationWithZoom(0);
		script += scriptGenerator.drawPolyhedron();
		script += scriptGenerator.drawAxes();
		script += "draw axes* on; draw poly* on; ";
		script += "save STATE state_1;";
		
		jmol.evalString(script);
		
	}
}
