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
 * Quick demo of how to call CE-Symm programmatically.
 *
 * @author spencer
 *
 */
public class AleixDemoMC {

	public static void main(String[] args){

	
		AtomCache cache = new AtomCache();

		String name = "4dou";
		
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75);


		CeSymm ceSymm = new CeSymm();

		try {
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);

			AFPChain afpChain = ceSymm.align(ca1, ca2);
			afpChain.setName1(name);
			afpChain.setName2(name);
			int symOrder = afpChain.getBlockNum();
			
			MultipleAFP multipleAFP = new MultipleAFP(afpChain,ca1,symOrder);
			afpChain = multipleAFP.getAfpChain();
			//Where a new AFPChain is initialized that the block colors are reset?
			Color[] colors = ColorBrewer.Set1.getColorPalette(symOrder);
			afpChain.setBlockColors(colors);
			
			SymmetryJmol jmol = SymmetryDisplay.display(afpChain, ca1, ca2);
						
			//int symmNr = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
			//System.out.println("Symmetry order of: " + symmNr);
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
