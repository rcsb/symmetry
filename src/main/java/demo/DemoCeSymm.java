package demo;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ReducedChemCompProvider;
import org.biojava.nbio.structure.scop.ScopFactory;

/**
 * Quick demo of how to call CE-Symm programmatically.
 *
 * @author spencer
 *
 */
public class DemoCeSymm {

	public static void main(String[] args){

	
		AtomCache cache = new AtomCache();

		String name = "4dou";
		name = "d1kcwa1";
		name = "3C1B.A:,B:,C:,D:,I:,E:,F:,G:,H:,J:"; // full nucleosome
		name = "3C1B.I:,J:"; // nucleosome DNA
		//name = "3C1B.A:,B:,E:,F:,C:,D:,G:,H:"; // nucleosome second axis
		name = "3F2Q";
		name = "d1wp5a_"; //special case C6
		name = "3EU9.A"; //hard cases non-closed: d1rmga_
		//name = "1N0R.A";  //Ankyrin: 1N0R.A, 3EHQ.A, 3EU9.A
		name = "1B3U.A";  //other repeats: 1xqr.A, 
		name = "d2ivza1";
		name = "1ITB.A";
		//name = "1vzw";
		//name = "d1flgb_";

		ScopFactory.setScopDatabase(ScopFactory.VERSION_2_0_4);


		CeSymm ceSymm = new CeSymm();
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		params.setRefineMethod(RefineMethod.MULTIPLE);
		//params.setMaxSymmOrder(3);
		params.setOptimization(true);
		params.setSymmetryType(SymmetryType.CLOSED);

		try {
			ChemCompGroupFactory.setChemCompProvider(new ReducedChemCompProvider());
			Structure struct = cache.getStructure(name);
			Atom[] ca1 = StructureTools.getRepresentativeAtomArray(struct);
			Atom[] ca2 = StructureTools.getRepresentativeAtomArray(struct.clone());
//			Atom[] ca2 = StructureTools.cloneAtomArray(ca1); // loses ligands
			AFPChain afpChain = ceSymm.align(ca1, ca2);
			afpChain.setName1(name);
			afpChain.setName2(name);
			
			System.out.println(AfpChainWriter.toDBSearchResult(afpChain));
			
			SymmetryJmol jmol = new SymmetryJmol(afpChain, ca1);
			
			RotationAxis axis = new RotationAxis(afpChain);
			jmol.evalString(axis.getJmolScript(ca1));
			
			int symmNr = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
			System.out.println("Symmetry order of: " + symmNr);
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
