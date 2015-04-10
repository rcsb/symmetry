package demo;
import java.awt.Color;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.gui.SymmetryDisplay;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.align.symm.subunit.MultipleAFP;
import org.jcolorbrewer.ColorBrewer;

/**
 * Demo for various scripting analysis.
 *
 * @author lafita
 *
 */
public class AleixDemo {

	public static void main(String[] args){

	
		AtomCache cache = new AtomCache();

		String name = "1vzw";
		
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75);
		
		try {
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);

			CeSymm ceSymm = new CeSymm();

			AFPChain afpChain = ceSymm.align(ca1, ca2);
			afpChain.setName1(name);
			afpChain.setName2(name);
			
			SymmetryJmol jmol = SymmetryDisplay.display(afpChain, ca1, ca2);
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
