import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.dbscan.FindMirrorSymmetries;
import org.biojava3.structure.utils.SymmetryTools;


public class ShowMirrorSymmetry {
	public static void main(String[] args){
		//String name1="1A25.B";
		//String name2="1A25.B";

		String name1="d1z7xw1";
		String name2="d1z7xw1";
		
		try {
			AtomCache cache = new AtomCache("/Users/andreas/WORK/PDB/",true);
			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);


			Atom[] ca2M = SymmetryTools.mirrorCoordinates(ca2);
			ca2M = SymmetryTools.duplicateMirrorCA2(ca2M);

			AFPChain afpChain = FindMirrorSymmetries.align(ca1,ca2M,name1,true);
			afpChain.setAlgorithmName(CeMain.algorithmName);

			//boolean isSignificant =  afp.isSignificantResult();
			boolean isSignificant = false;
			if ( afpChain.getTMScore() > 0.3 && afpChain.isSignificantResult()) {
				isSignificant = true;
				
			}
			
			StructureAlignmentDisplay.display(afpChain, ca1, ca2M);
			
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}

