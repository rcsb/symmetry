import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.PDBFileParser;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.dbscan.FindRotationSymmetries;
import org.biojava3.structure.utils.SymmetryTools;

/** Try to identify all possible symmetries by iterating resursively over all results and disabling the diagonal of each previous result.
 * 
 * @author andreas
 *
 */
public class IdentifyAllSymmetries {

	public static void main(String[] args){


		int fragmentLength = 8;
		AtomCache cache = new AtomCache("/Users/andreas/WORK/PDB/", true);

		// beta trefoil
		String name1 = "1JQZ.A";
		String name2 = "1JQZ.A";

		CeParameters params = new CeParameters();

		try {

			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2O = cache.getAtoms(name2);

			//ca2O = mirrorCoordinates(ca2O);

			Atom[] ca2 = SymmetryTools.duplicateCA2(ca2O);


			long timeS = System.currentTimeMillis();
			Matrix origM = null;
			
			AFPChain afpChain = new AFPChain();
			afpChain.setName1(name1);
			afpChain.setName2(name2);
			
			CECalculator calculator = new CECalculator(params);
			
			origM = align(afpChain,name1, name2, ca1,ca2, params, origM, calculator);
			//Matrix origM = dk.align(afpChain,chainId1, chainId2, ca1,ca2, params, null);

			long timeE = System.currentTimeMillis();
			System.out.println("===== Alignment took " + (timeE-timeS) + " ms.");
			StructureAlignmentDisplay.display(afpChain, ca1, ca2);


			double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
			afpChain.setTMScore(tmScore);
			System.out.println(AfpChainWriter.toScoresList(afpChain));
			afpChain.setDistanceMatrix((Matrix)origM.clone());
			double tmS = -1;
			double prevTm = -1;
			
			int i =0 ; ;
			//while (! String.format("%.4f",afpChain.getTMScore()).equals(String.format("%.4f", prevTm)) ){
			while ( ! String.format("%.4f",afpChain.getTMScore()).equals(String.format("%.4f", prevTm)) ){
				i++;
				System.out.println("Iteration : " + i + " " + afpChain.getTMScore());
				origM = align(afpChain,name1, name2, ca1,ca2, params, origM, calculator);
				AFPChain c = (AFPChain) afpChain.clone();
				StructureAlignmentDisplay.display(c, ca1, ca2);

				double tmScore2 = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
				afpChain.setTMScore(tmScore2);
				afpChain.setDistanceMatrix((Matrix)origM.clone());
				
				System.out.println(AfpChainWriter.toScoresList(afpChain));
				prevTm = tmS;
				tmS = tmScore2;
			}
			


		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public static Matrix align(AFPChain afpChain, String name1, String name2, Atom[] ca1, Atom[] ca2, 
			CeParameters params, Matrix origM, CECalculator calculator) throws StructureException{

		int fragmentLength = params.getWinSize();


		Atom[] ca2clone = SymmetryTools.cloneAtoms(ca2);



		////


		//System.out.println("rows "  + rows + " " + cols +
		//      " ca1 l " + ca1.length + " ca2 l " + ca2.length);

		//double[] dist1 = AlignTools.getDiagonalAtK(ca1, k);

		//double[] dist2 = AlignTools.getDiagonalAtK(ca2, k);

		int rows = ca1.length ;
		int cols = ca2.length ;

		// Matrix that tracks similarity of a fragment of length fragmentLength
		// starting a position i,j.

		//Matrix diffDistMax = getDkMatrix(ca1, ca2, fragmentLength, dist1, dist2, rows, cols);
		//showMatrix(diffDistMax,"Dk approach for setting initial alignment matrix");

		

		if ( origM == null) {
			//Build alignment ca1 to ca2-ca2

			afpChain.setName1(name1);
			afpChain.setName2(name2);
			afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);

			origM = SymmetryTools.blankOutCEOrig(ca2, rows, cols, calculator, origM);

			//SymmetryTools.showMatrix(origM, "original CE matrix");
			


		} else {
			
					
			// we are doing an iteration on a previous alignment
			// mask the previous alignment
			//afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);
			origM =  SymmetryTools.blankOutPreviousAlignment(afpChain, ca2, rows, cols, calculator,origM);
			System.out.println("BLANK OUT PREVIOUS");
			//SymmetryTools.showMatrix(origM, "iteration  matrix");

		}

		Matrix clone =(Matrix) origM.clone();

		calculator.setMatMatrix(clone.getArray());

		calculator.traceFragmentMatrix( afpChain,ca1, ca2clone);
		
		calculator.nextStep( afpChain,ca1, ca2clone);

		// afpChain = CeMain.filterDuplicateAFPs(afpChain, calculator, ca1, ca2clone);
		afpChain.setAlgorithmName("CE-symmetry");
		//afpChain.setVersion("0.0000001");


		return origM;



	}
}