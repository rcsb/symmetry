package demo;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.align.symm.CeSymm;


/** an example of CeSymm that is using sequence conservation as part of the alignment
 * 
 * @author Andreas Prlic
 *
 */
public class DemoCeSymmSeqConserv {

	public static void main(String[] args){
		
		String tmpDir = System.getProperty("java.io.tmpdir");
		AtomCache cache = new AtomCache(tmpDir,tmpDir);

		String name = "d1jlya1";

		CeSymm ceSymm = new CeSymm();

		try {
			
			Atom[] ca1 = cache.getAtoms(name); 
			Atom[] ca2 = cache.getAtoms(name);

			
			CeParameters params = (CeParameters) ceSymm.getParameters();
			if ( params == null) {
				params = new CeParameters();
				ceSymm.setParameters(params);
			}
				
			// here how to change the aa subst matrix, SDM is the default matrix
			//String matrixName = "PRLA000101";
			//SubstitutionMatrix<AminoAcidCompound> sdm = SubstitutionMatrixHelper.getMatrixFromAAINDEX(matrixName);			
			
			
			//SubstitutionMatrix<AminoAcidCompound> max = SubstitutionMatrixHelper.getBlosum85();
			//params.setSubstitutionMatrix(max);			
			
			
			params.setSeqWeight(2.0);
						
			AFPChain afpChain = ceSymm.align(ca1, ca2, params);
			afpChain.setName1(name);
			afpChain.setName2(name);

			System.out.println(AfpChainWriter.toDBSearchResult(afpChain));

			StructureAlignmentDisplay.display(afpChain, ca1, ca2);

		} catch (Exception e){
			e.printStackTrace();
		}
	}
}