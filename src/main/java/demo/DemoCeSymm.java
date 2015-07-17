package demo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.ChainSorter;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.gui.SymmetryDisplay;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Quick demo of how to call CE-Symm programmatically.
 * Some examples of different symmetry are proposed.
 *
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class DemoCeSymm {

	public static void main(String[] args) 
			throws IOException, StructureException {

		/* 
		 * Some examples:
		 * 
		 * CLOSED
		 * 2-fold: 1hiv.A, 
		 * 3-fold: 4i4q, 4dou
		 * 5-fold: 2jaj.A
		 * 6-fold: 1u6d
		 * 7-fold: 1jof.A
		 * 8-fold: 1vzw, d1i4na_
		 * 
		 * OPEN
		 * ankyrin: 1n0r.A, 3ehq.A
		 * leucine repeats: 2bnh.A
		 * helical: 1d0b.A
		 * 
		 * MULTIPLE AXES
		 * dihedral: 4hhb, 1vym
		 * hierarchical: 4gcr, 1ppr.O, 1hiv
		 * 
		 * - For more examples see the symmetry benchmark
		 */

		//Set the name of the protein structure to analyze
		String name = "1u6d";
		List<Atom[]> atoms = new ArrayList<Atom[]>();

		//Download the atoms and sort them sequentially by chains
		AtomCache cache = new AtomCache();
		Atom[] ca = ChainSorter.cyclicSorter(cache.getStructure(name));
		atoms.add(ca);

		CeSymm ceSymm = new CeSymm();

		//Choose some parameters
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setSymmetryType(SymmetryType.AUTO);
		params.setOptimization(true);

		//Run the alignment
		MultipleAlignment symmetry = ceSymm.align(atoms);

		//Display the results in jmol
		SymmetryDisplay.display(symmetry);
	}
}
