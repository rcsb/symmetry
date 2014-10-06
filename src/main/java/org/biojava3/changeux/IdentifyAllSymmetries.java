package org.biojava3.changeux;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.AbstractStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.ce.MatrixListener;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.utils.SymmetryTools;

/** Try to identify all possible symmetries by iterating resursively over all results and disabling the diagonal of each previous result.
 * 
 * @author andreas
 * @deprecated
 */
@Deprecated
public class IdentifyAllSymmetries  implements MatrixListener {

	public static final String algorithmName = "jCE-SYMM";

	public static final String version = "1.0";
	
	
	AFPChain afpChain;
	List<AFPChain> prevAligs = new ArrayList<AFPChain>();
	Atom[] ca1;
	Atom[] ca2;
	int rows;
	int cols;
	CECalculator calculator;
	CeParameters params;
	int loopCount ;
	int maxNrAlternatives = Integer.MAX_VALUE;
	boolean displayJmol = true;
	
	
	public static void main(String[] args){

		String name1 = args[0];
		String name2 = args[1];
		//String pdbPath = args[2];

		ProgressBar progress = ProgressBar.createAndShowGUI();
		progress.actionPerformed(null);
		progress.addStatus("detecting symmetry in  " + name1 + " vs. " + name2);
		
		UserConfiguration config = WebStartMain.getWebStartConfig();
		
		String pdbPath = config.getPdbFilePath();
		boolean isSplit = config.isSplit();
			
		int fragmentLength = 8;
		AtomCache cache = new AtomCache(pdbPath,pdbPath, isSplit);

		//1JQZ.A -  beta trefoil
		//1MSO - insulin
		// 2PHL.A - Phaseolin 

		//String name1 = name1;
		//String name2 = name2;

		IdentifyAllSymmetries me = new IdentifyAllSymmetries();
		me.setDisplayJmol(true);
		List<AFPChain> results = me.indentifyAllSymmetries(name1, name2, cache, fragmentLength,progress);
		
		int pos = -1;
		String legend = "# pos\tname1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\ttmScore\taligLen" ;
		progress.addStatus(legend);
		for (AFPChain afpChain : results){
			pos++;
			
			progress.addStatus("# " + pos + "\t" + toDBSearchResult(afpChain));
		}
		progress.done();
	}

	public void setMaxNrAlternatives(int max){
		maxNrAlternatives = max;
	}
	
	public int getMaxNrAlternatives(){
		return maxNrAlternatives;
	}
	
	
	
	public boolean isDisplayJmol() {
		return displayJmol;
	}

	public void setDisplayJmol(boolean displayJmol) {
		this.displayJmol = displayJmol;
	}

	public static String toDBSearchResult(AFPChain afpChain)
	{
		StringBuffer str = new StringBuffer();

		str.append(afpChain.getName1());
		str.append("\t");
		str.append(afpChain.getName2());
		str.append("\t");
		str.append(String.format("%.2f",afpChain.getAlignScore()));
		str.append("\t");     		
		str.append(String.format("%.2f",afpChain.getProbability()));		
		str.append("\t");
		str.append(String.format("%.2f",afpChain.getTotalRmsdOpt()));
		str.append("\t");
		str.append(afpChain.getCa1Length());
		str.append("\t");
		str.append(afpChain.getCa2Length());      
		str.append("\t");
		str.append(afpChain.getCoverage1());
		str.append("\t");
		str.append(afpChain.getCoverage2());
		str.append("\t");
		str.append(String.format("%.2f",afpChain.getTMScore()));
		str.append("\t");
		str.append(afpChain.getOptLength());
		return str.toString();
	}

	public List<AFPChain> indentifyAllSymmetries(String name1, String name2, AtomCache cache, int fragmentLength , ProgressBar progress){
		loopCount = 0;
		params = new CeParameters();

		//params.setMaxOptRMSD(3.5);
		params.setWinSize(fragmentLength);
		prevAligs = new ArrayList<AFPChain>();
		try {

			ca1 = cache.getAtoms(name1);
			Atom[] ca2O = cache.getAtoms(name2);

			//ca2O = mirrorCoordinates(ca2O);

			ca2 = StructureTools.duplicateCA2(ca2O);
			rows = ca1.length ;
			cols = ca2.length ;

			long timeS = System.currentTimeMillis();
			Matrix origM = null;

			AFPChain myAFP = new AFPChain();
			myAFP.setName1(name1);
			myAFP.setName2(name2);

			calculator = new CECalculator(params);
			calculator.addMatrixListener(this);

			int i =0 ; ;
			//while (! String.format("%.4f",afpChain.getTMScore()).equals(String.format("%.4f", prevTm)) ){
			//while ( ! String.format("%.4f",afpChain.getTMScore()).equals(String.format("%.4f", prevTm)) ){
			while ( (afpChain == null || isSignificant(myAFP) ) &&  i < maxNrAlternatives) {
				if ( progress == null) {
					//System.out.println("Iteration : " + i + " " + myAFP.getTMScore());
				}
				else {
					progress.addStatus("Analyzing alternative alignment " + i);
				}
				//this.afpChain = (AFPChain) myAFP.clone();
				afpChain = myAFP;
				if ( origM != null) {
					myAFP.setDistanceMatrix((Matrix)origM.clone());
				}
				origM = align(myAFP,name1, name2, ca1,ca2, params, origM, calculator, i);



				//				if ( String.format("%.4f",myAFP.getTMScore()).equals(String.format("%.4f", prevTm)) ) {
				//					System.out.println("seems I have converged...");
				//					break;
				//				}

				//System.out.println(AfpChainWriter.toScoresList(afpChain));
				double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
				myAFP.setTMScore(tmScore2);
				
				i++;
				if ( isSignificant(myAFP)) {
					prevAligs.add((AFPChain) myAFP.clone());

				}
			}


			for ( AFPChain a : prevAligs){
				double tmScore2 = AFPChainScorer.getTMScore(a, ca1, ca2);
				a.setTMScore(tmScore2);
				
				a = CeCPMain.filterDuplicateAFPs(a, calculator, ca1, ca2);
				if ( displayJmol)
					showCurrentAlig(a, ca1, ca2);
			}

			if ( prevAligs.size()==0){
				double tmScore2 = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
				afpChain.setTMScore(tmScore2);
				try {
					afpChain = CeCPMain.filterDuplicateAFPs(afpChain, calculator, ca1, ca2);
				} catch (Exception e){
					return prevAligs;
				}
				
				//afpChain = CeMain.filterDuplicateAFPs(afpChain,calculator,ca1,ca2);
				prevAligs.add(afpChain);
				if ( displayJmol) {
					//if ( afpChain.getProbability() > 3.5) 
					showCurrentAlig(afpChain, ca1, ca2);
				}
			}
			//System.out.println("We went through " + i + " iterations. Final TM score: " + afpChain.getTMScore());
			//afpChain.setAlgorithmName("CE-symmetry final result ");
			//StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			//jmol.evalString("draw l1 line 100 {0 0 0} (1:A.CA/1) ; draw l2 line 100 {0 0 0} (1:A.CA/2);" );
		} catch (Exception e){
			e.printStackTrace();
		}
		return prevAligs;
	}

	public static  boolean isSignificant(AFPChain myAFP){

		//|| 
		return ( (myAFP.getTMScore() >= 0.35  && myAFP.getProbability() >= 3.5 ) && myAFP.getTotalRmsdOpt() < 5.0);
	}

	private void showCurrentAlig(AFPChain myAFP, Atom[] ca12, Atom[] ca22) throws StructureException {
		AFPChain c = (AFPChain) myAFP.clone();
		StructureAlignmentJmol jmol =  StructureAlignmentDisplay.display(c, ca1, ca2);

		
		// draw a line from center of gravity to N terminus
		
		ResidueNumber res1 = ca1[0].getGroup().getResidueNumber();
		ResidueNumber res2 = ca2[0].getGroup().getResidueNumber();
		String chainId1 = ca1[0].getGroup().getChain().getChainID();
		String chainId2 = ca2[0].getGroup().getChain().getChainID();
		
		Atom centroid1 = Calc.getCentroid(ca12);
		Atom centroid2 = Calc.getCentroid(ca22);
		
		String cs1 = "{"+centroid1.getX() + " " + centroid1.getY() + " " +centroid1.getZ()+"}";
		String cs2 = "{"+centroid2.getX() + " " + centroid2.getY() + " " +centroid2.getZ()+"}";
		
		jmol.evalString("draw l1 line 100 "+cs1+" (" + 
				res1.getSeqNum()+":" + chainId1+".CA/1) ; draw l2 line 100 "+cs2+" ("+
				res2.getSeqNum()+":" + chainId2+".CA/2);" );


	}



	private static Matrix align(AFPChain afpChain, String name1, String name2, Atom[] ca1, Atom[] ca2, 
			CeParameters params, Matrix origM, CECalculator calculator, int counter) throws StructureException{

		int fragmentLength = params.getWinSize();
		//if ( ca1.length > 200 && ca2.length > 200 )
		//	fragmentLength = 25;

		//params.setWinSize(fragmentLength);

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


		int blankWindowSize = fragmentLength ;
		if ( origM == null) {
			//Build alignment ca1 to ca2-ca2

			afpChain.setName1(name1);
			afpChain.setName2(name2);
			afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);

			origM =  SymmetryTools.blankOutPreviousAlignment(afpChain, ca2, rows, cols, calculator, null, blankWindowSize);

			//SymmetryTools.showMatrix(origM, "original CE matrix");



		} else {

			// we are doing an iteration on a previous alignment
			// mask the previous alignment
			//afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);
			origM =  SymmetryTools.blankOutPreviousAlignment(afpChain, ca2, rows, cols, calculator, origM, blankWindowSize);
			//afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);
			//System.out.println("BLANK OUT PREVIOUS");
			//SymmetryTools.showMatrix(origM, "iteration  matrix " +counter);

		}		

		Matrix clone =(Matrix) origM.clone();

		// that's the matrix to run the alignment on..
		calculator.setMatMatrix(clone.getArray());

		calculator.traceFragmentMatrix( afpChain,ca1, ca2clone);

		calculator.nextStep( afpChain,ca1, ca2clone);

		// afpChain = CeMain.filterDuplicateAFPs(afpChain, calculator, ca1, ca2clone);
		afpChain.setAlgorithmName("CE-symmetry step " + counter);

		//afpChain.setVersion("0.0000001");

		//if ( isSignificant(afpChain)) {
		//AFPChain c = (AFPChain) afpChain.clone();
		//afpChain = CeMain.filterDuplicateAFPs(c,calculator,ca1,ca2);
		//afpChain.setDistanceMatrix(origM);
		//	afpChain = CeMain.filterDuplicateAFPs(afpChain, calculator, ca1, ca2clone);
		//afpChain.setDistanceMatrix((Matrix)origM.clone());
		//	}

		return origM;

	}

	public double[][] matrixInOptimizer(double[][] max) {
		//System.out.println("In optimizer masking known alignments " );

		int fragmentLength = params.getWinSize();
		//if ( ca1.length > 100 && ca2.length > 100 )
		//	fragmentLength = fragmentLength * 2;
		int blankWindowSize = fragmentLength ;

		Matrix origM = new Matrix(max);
		//SymmetryTools.showMatrix((Matrix)origM.clone(), "before mask "  + loopCount  );
		//SymmetryTools.showMatrix(origM, "iteration  matrix " + loopCount + " before");

		//System.out.println("iteration X..." + loopCount);

		for ( AFPChain prevAlig : prevAligs) {

			origM =  SymmetryTools.blankOutPreviousAlignment(prevAlig, 
					ca2, rows, cols, calculator, origM, blankWindowSize);
		}

		if ( afpChain != null) {
			// blank out current alignment
			origM = SymmetryTools.blankOutPreviousAlignment(afpChain, 
					ca2, rows, cols, calculator, origM, blankWindowSize);
		}

		//origM =  SymmetryTools.blankOutPreviousAlignment(afpChain, 
		//		   ca2, rows, cols, calculator, origM, blankWindowSize);

		max = origM.getArray();

		//SymmetryTools.showMatrix((Matrix)origM.clone(), "in optimizer "  + loopCount  );
		//SymmetryTools.showMatrix(origM, "iteration  matrix " + loopCount + " after");
		loopCount++;
		return max;
	}



	public boolean[][] initializeBreakFlag(boolean[][] breakFlag) {
		int fragmentLength = params.getWinSize();
		try {
			if ( afpChain != null) {
				//afpChain = CeMain.filterDuplicateAFPs(afpChain, calculator, ca1, ca2);
				breakFlag = SymmetryTools.blankOutBreakFlag(afpChain, 
						ca2, rows, cols, calculator, breakFlag, fragmentLength);
			}

			for ( AFPChain prevAlig : prevAligs) {
				prevAlig = CeCPMain.filterDuplicateAFPs(prevAlig, calculator, ca1, ca2);
				breakFlag =  SymmetryTools.blankOutBreakFlag(prevAlig, 
						ca2, rows, cols, calculator, breakFlag, fragmentLength);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		return breakFlag;

	}

	public CECalculator getCalculator() {
		return calculator;
	}

	public void setCalculator(CECalculator calculator) {
		this.calculator = calculator;
	}

	

	
	
	
}