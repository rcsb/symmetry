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
			
			String chainId1 = "1AUO.A";
			String chainId2 = "1AUO.A";


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

	
	

	
	

	
	
}