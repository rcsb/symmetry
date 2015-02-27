package demo;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;

/**
 * Demo for the CE-symm with recursive symmetry detection.
 * Given the order of symmetry 
 * 
 * @author aleix
 *
 */
public class AleixDemoCeSymm {

	public static void main(String[] args){

	
		AtomCache cache = new AtomCache();

		String name = "4HHB";
		
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75);

		CeSymm ceSymm = new CeSymm();

		try {
			
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);
			
			CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
			params.setMaxNrAlternatives(2);
			

			AFPChain afpChain = ceSymm.align(ca1, ca2, params);
			afpChain.setName1(name);
			afpChain.setName2(name);
				
			System.out.println(AfpChainWriter.toDBSearchResult(afpChain));
			
			StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			
			RotationAxis axis = new RotationAxis(afpChain);
			jmol.evalString(axis.getJmolScript(ca1));
			
			int symmNr = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
			System.out.println("Symmetry order of: " + symmNr);
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
