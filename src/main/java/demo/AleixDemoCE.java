// 24 February 2015 -- Aleix Lafita
// Biojava 3 Tutorial on Github - Structural Alignment

package demo;

import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.ce.GuiWrapper;
import org.biojava.nbio.structure.align.gui.AlignmentGui;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;


public class AleixDemoCE {

	public static void main(String[] args) {
		
		try{
		//Fetch CA atoms for the structures to be aligned
		String name1 = "3cna.A";
		String name2 = "2pel";
		AtomCache cache = new AtomCache();
		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);
		
		//Get StructureAlignment instance
		StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
		
		//Perform the alignment
		AFPChain afpChain = algorithm.align(ca1,ca2);
		
		//Print text output
		System.out.println(afpChain.toCE(ca1,ca2));
		
		//Display the alignment in Jmol
		GuiWrapper.display(afpChain,ca1,ca2);
		
		} catch (Exception e){
			e.printStackTrace();
		}
				
	}

}
