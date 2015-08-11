package demo;

import java.io.IOException;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.align.symm.CeSymm;
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
		 * - For more examples see the symmetry benchmark
		 */

		//Set the name of the protein structure to analyze
		String name = "1hiv.A";

		//Download the atoms and sort them sequentially by chains
		AtomCache cache = new AtomCache();
		Structure s = cache.getStructure(name);
		Atom[] ca1 = StructureTools.getRepresentativeAtomArray(s);
		Atom[] ca2 = StructureTools.getRepresentativeAtomArray(s.clone());

		CeSymm ceSymm = new CeSymm();

		//Choose some parameters
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setSymmetryType(SymmetryType.AUTO);

		//Run the alignment
		AFPChain symmetry = ceSymm.align(ca1, ca2);

		//Display the results in jmol
		StructureAlignmentDisplay.display(symmetry, ca1, ca2);
	}
}
