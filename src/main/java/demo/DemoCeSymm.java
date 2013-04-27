package demo;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.SymmRefiner;


public class DemoCeSymm {

	public static void main(String[] args){

	
		AtomCache cache = new AtomCache();

		String name = "d1jlya1";


		CeSymm ceSymm = new CeSymm();

		try {
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);

			AFPChain afpChain = ceSymm.align(ca1, ca2);
			afpChain.setName1(name);
			afpChain.setName2(name);
			
			System.out.println(AfpChainWriter.toDBSearchResult(afpChain));
			
			StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			
			RotationAxis axis = new RotationAxis(afpChain);
			jmol.evalString(axis.getJmolScript(ca1));
			
			int symmNr = CeSymm.getSymmetryOrder(afpChain);
			System.out.println("Symmetry order of: " + symmNr);
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
