import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.dbscan.FindMirrorSymmetries;
import org.biojava3.structure.utils.SymmetryTools;


public class ShowMirrorSymmetry {
	public static void main(String[] args){
		//String name1="1A25.B";
		//String name2="1A25.B";

		String name1="3QGM.D";
		String name2="3QGM.D";
		
		try {
			AtomCache cache = new AtomCache("/Users/andreas/WORK/PDB/",true);
			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);


			Atom[] ca2M = SymmetryTools.mirrorCoordinates(ca2);
			ca2M = SymmetryTools.duplicateMirrorCA2(ca2M);

			AFPChain afp = FindMirrorSymmetries.align(ca1,ca2M,name1);
			afp.setAlgorithmName(CeMain.algorithmName);

			//boolean isSignificant =  afp.isSignificantResult();
			boolean isSignificant = false;
			if ( afp.getTMScore() > 0.3 && afp.isSignificantResult()) {
				isSignificant = true;
				
			}
			
			StructureAlignmentDisplay.display(afp, ca1, ca2M);
			
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}

