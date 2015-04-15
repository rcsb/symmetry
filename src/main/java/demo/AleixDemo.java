package demo;
import java.io.IOException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.gui.SymmetryDisplay;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;

/**
 * Demo for the CeSymm alignment and display.
 *
 * @author lafita
 * 
 */
public class AleixDemo {

	public static void main(String[] args) throws StructureException, IOException{

	
		AtomCache cache = new AtomCache();
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
		//ScopFactory.setScopDatabase(ScopFactory.VERSION_2_0_1, true);
		
		//Easy cases: 4i4q, 4dou
		//Hard cases: d2vdka_,d1n6dd3
		//Better MULTIPLE: 2i5i.a
		String name = "d1n6dd3";
	
		Atom[] ca1 = cache.getAtoms(name);
		Atom[] ca2 = cache.getAtoms(name);
		System.out.println(ca1.length);

		CeSymm ceSymm = new CeSymm();
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		params.setRefineMethod(RefineMethod.MONTE_CARLO);
		AFPChain afpChain = ceSymm.align(ca1, ca2);
		
		afpChain.setName1(name);
		afpChain.setName2(name);
		
		SymmetryJmol jmol = SymmetryDisplay.display(afpChain, ca1, ca2);
			
	}
}
