import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.align.symm.CeSymm;


public class DemoCeSymm {

	public static void main(String[] args){

		String tmpDir = System.getProperty("java.io.tmpdir");
		AtomCache cache = new AtomCache(tmpDir, true);

		String name = "d1jlya1";

		CeSymm ceSymm = new CeSymm();

		try {
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);

			AFPChain afpChain = ceSymm.align(ca1, ca2);
			afpChain.setName1(name);
			afpChain.setName2(name);
			
			System.out.println(AfpChainWriter.toDBSearchResult(afpChain));
			
			StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
