import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JFrame;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.helper.AlignTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.gui.ScaleableMatrixPanel;
import org.biojava.bio.structure.io.PDBFileParser;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.dbscan.FindRotationSymmetries;
import org.biojava3.structure.utils.SymmetryTools;



public class ShowRotationSymmetry {


	public static void main(String[] args){
		


		int fragmentLength = 8;


		try {
			AtomCache cache = new AtomCache("/Users/andreas/WORK/PDB/", true);
			Logger.getLogger(PDBFileParser.class.getName()).setLevel(Level.INFO);

			// big one
			//String chainId1 = "1A9X.G";
			//String chainId2 = "1A9X.G";


			//String chainId1 = "3QGM.D";
			//String chainId2 = "3QGM.D";

			//String chainId1 = "1A6S.A";
			//String chainId2 = "1A6S.A";
			
			String chainId1 = "d1z7xw1";
			String chainId2 = "d1z7xw1";


			// intra and intermolecularsymm
			//String chainId1 = "1jnr.D";
			//String chainId2 = "1jnr.B";

			// intramolecularsymm
			//String chainId1 = "1mer.A";
			//String chainId2 = "1mer.A";


			// beta trefoil
			//String chainId1 = "1JQZ.A";
			//String chainId2 = "1JQZ.A";

			// four helix bundle
			//String chainId1 = "2MHR.A";
			//String chainId2 = "2MHR.A";

			// symm
			//String chainId1 = "1bp7.A";
			//String chainId2 = "1bp7.B";


			//String chainId1 = "1BYI.A";
			//String chainId2 = "1BYI.A";


			//String chainId1 = "1FE0.A";
			//String chainId2 = "1FE0.A";

			//String chainId1 = "1BTL.A";
			//String chainId2 = "1BTL.A";

			//Structure s = cache.getStructure(chainId1);
			//System.out.println("SITES: " + s.getSites());
			//StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			//jmol.setStructure(s);

			Atom[] ca1 = cache.getAtoms(chainId1);
			Atom[] ca2O = cache.getAtoms(chainId2);

			//ca2O = mirrorCoordinates(ca2O);

			Atom[] ca2 = SymmetryTools.duplicateCA2(ca2O);

			CeParameters params = new CeParameters();
			params.setWinSize(fragmentLength);
			// set the maximum gap size to unlimited 



			long timeS = System.currentTimeMillis();
			
			AFPChain afpChain = FindRotationSymmetries.align(ca1, ca2, chainId1, true);
			//Matrix origM = dk.align(afpChain,chainId1, chainId2, ca1,ca2, params, null);
			
			long timeE = System.currentTimeMillis();
			System.out.println("===== Alignment took " + (timeE-timeS) + " ms.");
			StructureAlignmentDisplay.display(afpChain, ca1, ca2);


			double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
			afpChain.setTMScore(tmScore);
			System.out.println(AfpChainWriter.toScoresList(afpChain));

//			int i =0 ; ;
//			while ( afpChain.getTMScore() > 0.3){
//				i++;
//				System.out.println("Iteration : " + i + " " + afpChain.getTMScore());
//				origM = dk.align(afpChain,chainId1, chainId2, ca1,ca2, params, origM);
//				StructureAlignmentDisplay.display(afpChain, ca1, ca2);
//				double tmScore2 = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
//				afpChain.setTMScore(tmScore2);
//			}
			//long timeS2 = System.currentTimeMillis();
			//AFPChain afpChain2 = new CeMain().align(ca1,ca2,params);            
			//long timeE2 = System.currentTimeMillis();
			//System.out.println("===== Alignment took " + (timeE2-timeS2) + " ms.");



			//StructureAlignmentDisplay.display(afpChain2, ca1, ca2O);


			//double tmScore2 = AFPChainScorer.getTMScore(afpChain2, ca1, ca2);
			//afpChain2.setTMScore(tmScore2);
			//System.out.println(AfpChainWriter.toScoresList(afpChain2));


		} catch (Exception e){
			e.printStackTrace();
		}
	}

	
	

	
	public Matrix align(AFPChain afpChain, String name1, String name2, Atom[] ca1, Atom[] ca2, 
			CeParameters params, Matrix origM) throws StructureException{


		int k = 6;
		int fragmentLength = params.getWinSize();


		Atom[] ca2clone = SymmetryTools.cloneAtoms(ca2);



		////


		//System.out.println("rows "  + rows + " " + cols +
		//      " ca1 l " + ca1.length + " ca2 l " + ca2.length);

		//double[] dist1 = AlignTools.getDiagonalAtK(ca1, k);

		//double[] dist2 = AlignTools.getDiagonalAtK(ca2, k);

		int rows = ca1.length - fragmentLength - k + 1;
		int cols = ca2.length - fragmentLength - k + 1;

		// Matrix that tracks similarity of a fragment of length fragmentLength
		// starting a position i,j.

		//Matrix diffDistMax = getDkMatrix(ca1, ca2, fragmentLength, dist1, dist2, rows, cols);
		//showMatrix(diffDistMax,"Dk approach for setting initial alignment matrix");

		CECalculator calculator = new CECalculator(params);

		if ( origM == null) {
			//Build alignment ca1 to ca2-ca2

			afpChain.setName1(name1);
			afpChain.setName2(name2);
			afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);

			origM = SymmetryTools.blankOutCEOrig(ca2, rows, cols, calculator, origM);

			SymmetryTools.showMatrix(origM, "original CE matrix");



		} else {
			// we are doing an iteration on a previous alignment
			// mask the previous alignment
			afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);
			origM =  SymmetryTools.blankOutPreviousAlignment(afpChain, ca2, rows, cols, calculator,origM);

			SymmetryTools.showMatrix(origM, "iteration  matrix");

		}

		Matrix clone =(Matrix) origM.clone();
		calculator.setMatMatrix(clone.getArray());

		//calculator.setMatMatrix(diffDistMax.getArray());

		//afpChain.setDistanceMatrix();

		calculator.traceFragmentMatrix( afpChain,ca1, ca2clone);
		calculator.nextStep( afpChain,ca1, ca2clone);

		// afpChain = CeMain.filterDuplicateAFPs(afpChain, calculator, ca1, ca2clone);
		afpChain.setAlgorithmName("CE-symmetry");
		//afpChain.setVersion("0.0000001");


		return origM;



	}

	
	
}