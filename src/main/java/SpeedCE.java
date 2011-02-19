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



public class SpeedCE {


	public static void main(String[] args){
		SpeedCE dk = new SpeedCE();


		int fragmentLength = 8;


		try {
			AtomCache cache = new AtomCache("/tmp/", true);
			Logger.getLogger(PDBFileParser.class.getName()).setLevel(Level.INFO);

			// big one
			//String chainId1 = "1A9X.G";
			//String chainId2 = "1A9X.G";


			String chainId1 = "3QGM.D";
			String chainId2 = "3QGM.D";


			//String chainId1 = "1ubi";
			//String chainId2 = "1ubi";


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

			Atom[] ca2 = duplicateCA2(ca2O);

			CeParameters params = new CeParameters();
			params.setWinSize(fragmentLength);
			// set the maximum gap size to unlimited 



			long timeS = System.currentTimeMillis();
			AFPChain afpChain = new AFPChain();
			Matrix origM = dk.align(afpChain,chainId1, chainId2, ca1,ca2, params, null);
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

	private static Atom[] mirrorCoordinates(Atom[] ca2O) {
		for(int i=0;i<ca2O.length;i++) {
			//ca2O[i].setX(-ca2O[i].getX());
			Group g = ca2O[i].getGroup();
			for ( Atom a : g.getAtoms()){
				a.setX(-a.getX());
			}
		}

		return ca2O;
	}

	public static Atom[] duplicateMirrorCA2(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length*2];

		int pos = ca2clone.length - 1;

		Chain c = new ChainImpl();
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos--;
		}


		// Duplicate ca2!
		for (Atom a : ca2){
			Group g = (Group)a.getGroup().clone();
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos--;
		}

		return ca2clone;


	}

	public static Atom[] duplicateCA2(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length*2];

		int pos = 0;

		Chain c = new ChainImpl();
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos++;
		}


		// Duplicate ca2!
		for (Atom a : ca2){
			Group g = (Group)a.getGroup().clone();
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos++;
		}

		return ca2clone;


	}

	public Matrix align(AFPChain afpChain, String name1, String name2, Atom[] ca1, Atom[] ca2, 
			CeParameters params, Matrix origM) throws StructureException{


		int k = 6;
		int fragmentLength = params.getWinSize();


		///
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length];

		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom

			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos++;
		}



		////


		//System.out.println("rows "  + rows + " " + cols +
		//      " ca1 l " + ca1.length + " ca2 l " + ca2.length);

		double[] dist1 = AlignTools.getDiagonalAtK(ca1, k);

		double[] dist2 = AlignTools.getDiagonalAtK(ca2, k);

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

			origM = blankOutCEOrig(ca2, rows, cols, calculator, origM);

			showMatrix(origM, "original CE matrix");



		} else {
			// we are doing an iteration on a previous alignment
			// mask the previous alignment
			afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);
			origM =  blankOutPreviousAlignment(afpChain, ca2, rows, cols, calculator,origM);

			showMatrix(origM, "iteration  matrix");

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

	private Matrix blankOutPreviousAlignment(AFPChain afpChain, Atom[] ca2,
			int rows, int cols, CECalculator calculator, Matrix max) {

		
		double[][] dist1 = calculator.getDist1();
		double[][] dist2 = calculator.getDist2();
		
		int[][][] optAln = afpChain.getOptAln();
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();
		for(int bk = 0; bk < blockNum; bk ++)       {

			//Matrix m= afpChain.getBlockRotationMatrix()[bk];
			//Atom shift = afpChain.getBlockShiftVector()[bk];
			for ( int i=0;i< optLen[bk];i++){
				int pos1 = optAln[bk][0][i];
				int pos2 = optAln[bk][1][i];
				// blank out area around these positions...

				int dist = 10 ;
				int start1 = Math.max(pos1-dist,0);
				int start2 = Math.max(pos2-dist,0);
				int end1 = Math.min(pos1+dist, rows-1);
				int end2 = Math.min(pos2+dist, cols-1);
				//System.out.println(pos1 + "  " + pos2 + " " + start1 + " " + end1 + " " + start2 + " " + end2);
				for ( int i1 = start1; i1< end1 ; i1++){
					dist1[i1][i1] = 99;
					
					for ( int j2 = start2 ; j2 < end2 ; j2++){
						//System.out.println(i1 + " " + j2);
						max.set(i1,j2,99);
						dist2[j2][j2] = 99;
					}
				}

			}

		}
		return max;

	}

	private Matrix blankOutCEOrig(Atom[] ca2, int rows, int cols,
			CECalculator calculator, Matrix origM) {
		origM = new Matrix( calculator.getMatMatrix());
		// symmetry hack, disable main diagonale

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				if ( diff < 15 ){
					origM.set(i,j, 99);
				}
				int diff2 = Math.abs(i-(j-ca2.length/2));
				if ( diff2 < 15 ){
					origM.set(i,j, 99);
				}
			}
		}
		return origM;
	}

	private Matrix  getDkMatrix(Atom[] ca1, Atom[] ca2,int fragmentLength,
			double[] dist1, double[] dist2, int rows, int cols) {
		Matrix diffDistMax =  Matrix.identity(ca1.length, ca2.length);

		for ( int i = 0 ; i< rows; i++){
			double score1 = 0;
			for ( int x=0 ; x < fragmentLength ; x++){
				score1 += dist1[i+x];
			}
			for ( int j = 0 ; j < cols ; j++){
				double score2 = 0;
				for ( int y=0 ; y < fragmentLength ; y++){
					score2 += dist2[j+y];
				}

				// if the intramolecular distances are very similar
				// the two scores should be similar, i.e. the difference is close to 0
				diffDistMax.set(i,j, Math.abs(score1-score2));
			}
		}


		// symmetry hack, disable main diagonale

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				if ( diff < 15 ){
					diffDistMax.set(i,j, 99);
				}
				int diff2 = Math.abs(i-(j-ca2.length/2));
				if ( diff2 < 15 ){
					diffDistMax.set(i,j, 99);
				}
			}
		}
		return diffDistMax;

	}

	private void showMatrix(Matrix m, String string) {
		ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
		JFrame frame = new JFrame();
		frame.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				JFrame f = (JFrame) e.getSource();
				f.setVisible(false);
				f.dispose();
				System.exit(0);
			}
		});

		smp.setMatrix((Matrix)m.clone());
		//smp.getMatrixPanel().setScale(0.8f);

		frame.setTitle(string);
		frame.getContentPane().add(smp);

		frame.pack();
		frame.setVisible(true);

	}
}