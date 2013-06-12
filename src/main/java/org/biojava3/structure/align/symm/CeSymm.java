package org.biojava3.structure.align.symm;
import java.io.IOException;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.AbstractStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.ce.MatrixListener;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AlignmentTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.utils.SymmetryTools;


/** Try to identify all possible symmetries by iterating resursively over all results and disabling the diagonal of each previous result.
 * 
 * @author andreas
 *
 */
public class CeSymm extends AbstractStructureAlignment implements MatrixListener, StructureAlignment {
	static final boolean debug = false;
	boolean displayJmol = false;

	public static final String algorithmName = "jCE-symmetry";

	public static final String version = "1.0";

	/**
	 * The penalty x residues from the main diagonal is: (standard score) - GRADIENT_EXP_COEFF*e^(-x) - GRADIENT_POLY_COEFF[n]*x^(-n) - GRADIENT_POLY_COEFF[n-1]*x^(-n+1) - GRADIENT_POLY_COEFF[n-2]*x^(-n+2) - ... - GRADIENT_POLY_COEFF[0]*x^0
	 */
	public static double[] GRADIENT_POLY_COEFF = {Integer.MIN_VALUE}; // from left to right: ..., quintic, quartic, cubic, quadratic, linear, constant; can be any length
	public static double GRADIENT_EXP_COEFF = 0; // the corresponding radix is e

	AFPChain afpChain;

	Atom[] ca1;
	Atom[] ca2;
	int rows;
	int cols;
	CECalculator calculator;
	CeParameters params;
	//int loopCount ;
	int maxNrAlternatives = 1;



	public static void main(String[] args){

		// used only for printing help...
		CeSymm ce = new CeSymm();

		if ( args.length < 2 ) {

			System.out.println(ce.printHelp());
			return;
		}

		if (args.length  == 0 ) {			
			System.out.println(ce.printHelp());
			return;			
		}

		if ( args.length == 1){
			if (args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("-help")|| args[0].equalsIgnoreCase("--help")){
				System.out.println(ce.printHelp());								
				return;
			}

		}

		CeSymmUserArgumentProcessor processor = new CeSymmUserArgumentProcessor(); //Responsible for creating a CeMain instance

		processor.process(args);

		//		String name1 = args[0];
		//		String name2 = args[1];
		//		//String pdbPath = args[2];
		//
		//		ProgressBar progress = ProgressBar.createAndShowGUI();
		//		progress.actionPerformed(null);
		//		progress.addStatus("detecting symmetry in  " + name1 + " vs. " + name2);
		//
		//		UserConfiguration config = WebStartMain.getWebStartConfig();
		//
		//		String pdbPath = config.getPdbFilePath();
		//		boolean isSplit = config.isSplit();
		//
		//		int fragmentLength = 8;
		//		AtomCache cache = new AtomCache(pdbPath, isSplit);
		//
		//		//1JQZ.A -  beta trefoil
		//		//1MSO - insulin
		//		// 2PHL.A - Phaseolin 
		//
		//		//String name1 = name1;
		//		//String name2 = name2;
		//
		//		CeSymm me = new CeSymm();
		//		List<AFPChain> results = new ArrayList<AFPChain>();
		//		try {
		//			results = me.indentifyAllSymmetries(name1, name2, cache, fragmentLength,progress);
		//		} catch (Exception e){
		//			e.printStackTrace();
		//		}
		//		int pos = -1;
		//		String legend = "# pos\tname1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\ttmScore\taligLen" ;
		//		progress.addStatus(legend);
		//		for (AFPChain afpChain : results){
		//			pos++;
		//
		//			progress.addStatus("# " + pos + "\t" + toDBSearchResult(afpChain));
		//		}
		//		progress.done();
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

	public AFPChain indentifyAllSymmetries(String name1, String name2, AtomCache cache, int fragmentLength ) throws StructureException,IOException{
		//loopCount = 0;
		params = new CeParameters();

		//params.setMaxOptRMSD(3.5);
		params.setWinSize(fragmentLength);


		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);


		return align(ca1,ca2,params);


	}

	public static  boolean isSignificant(AFPChain myAFP){

		//|| 
		return ( (myAFP.getTMScore() >= 0.35  || myAFP.getProbability() >= 3.5 ) && myAFP.getTotalRmsdOpt() < 5.0);
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



	private static Matrix align(AFPChain afpChain,  Atom[] ca1, Atom[] ca2, 
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

			afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);

			origM =  SymmetryTools.grayOutPreviousAlignment(afpChain, ca2, rows, cols, calculator, null, blankWindowSize, GRADIENT_POLY_COEFF, GRADIENT_EXP_COEFF);

			//SymmetryTools.showMatrix(origM, "original CE matrix");



		} else {

			// we are doing an iteration on a previous alignment
			// mask the previous alignment
			//afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);
			origM =  SymmetryTools.grayOutPreviousAlignment(afpChain, ca2, rows, cols, calculator, origM, blankWindowSize, GRADIENT_POLY_COEFF, GRADIENT_EXP_COEFF);
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
		//afpChain.setAlgorithmName("CE-symmetry step " + counter);
		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);

		//afpChain.setVersion("0.0000001");

		//if ( isSignificant(afpChain)) {
		//AFPChain c = (AFPChain) afpChain.clone();
		//afpChain = CeMain.filterDuplicateAFPs(c,calculator,ca1,ca2);
		afpChain.setDistanceMatrix(origM);
		//	afpChain = CeMain.filterDuplicateAFPs(afpChain, calculator, ca1, ca2clone);
		//afpChain.setDistanceMatrix((Matrix)origM.clone());
		//	}

		return origM;

	}

	public double[][] matrixInOptimizer(double[][] max) {

		return CECalculator.updateMatrixWithSequenceConservation(max, ca1, ca2, params);
	}



	public boolean[][] initializeBreakFlag(boolean[][] breakFlag) {
		int fragmentLength = params.getWinSize();
		try {
			if ( afpChain != null) {
				//afpChain = CeMain.filterDuplicateAFPs(afpChain, calculator, ca1, ca2);
				breakFlag = SymmetryTools.blankOutBreakFlag(afpChain, 
						ca2, rows, cols, calculator, breakFlag, fragmentLength);
			}

			//			for ( AFPChain prevAlig : prevAligs) {
			//				prevAlig = CeMain.filterDuplicateAFPs(prevAlig, calculator, ca1, ca2);
			//				breakFlag =  SymmetryTools.blankOutBreakFlag(prevAlig, 
			//						ca2, rows, cols, calculator, breakFlag, fragmentLength);
			//			}
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

	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {

		if (params == null)
			params = new CeParameters();

		return align(ca1,ca2,params);
	}

	public AFPChain align(Atom[] ca1, Atom[] ca2O, Object param)
			throws StructureException {
		if ( ! (param instanceof CeParameters))
			throw new IllegalArgumentException("CE algorithm needs an object of call CeParameters as argument.");

		this.params = (CeParameters) param;

		this.ca1=ca1;
		this.ca2=ca2O;

		//prevAligs = new ArrayList<AFPChain>();
		try {

			//ca2O = mirrorCoordinates(ca2O);

			ca2 = StructureTools.duplicateCA2(ca2O);
			rows = ca1.length ;
			cols = ca2.length ;

			Matrix origM = null;

			AFPChain myAFP = new AFPChain();

			calculator = new CECalculator(params);
			calculator.addMatrixListener(this);

			int i =0 ; ;
			//while (! String.format("%.4f",afpChain.getTMScore()).equals(String.format("%.4f", prevTm)) ){
			//while ( ! String.format("%.4f",afpChain.getTMScore()).equals(String.format("%.4f", prevTm)) ){
			while ( (afpChain == null || isSignificant(myAFP) ) &&  i < maxNrAlternatives) {

				//this.afpChain = (AFPChain) myAFP.clone();
				afpChain = myAFP;
				if ( origM != null) {
					myAFP.setDistanceMatrix((Matrix)origM.clone());
				}
				origM = align(myAFP, ca1,ca2, params, origM, calculator, i);



				//				if ( String.format("%.4f",myAFP.getTMScore()).equals(String.format("%.4f", prevTm)) ) {
				//					System.out.println("seems I have converged...");
				//					break;
				//				}

				//System.out.println(AfpChainWriter.toScoresList(afpChain));
				double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
				myAFP.setTMScore(tmScore2);

				i++;
				//				if ( isSignificant(myAFP)) {
				//					prevAligs.add((AFPChain) myAFP.clone());
				//
				//				}
			}



			double tmScore2 = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
			afpChain.setTMScore(tmScore2);
			try {
				afpChain = CeCPMain.postProcessAlignment(afpChain, ca1, ca2, calculator);
				afpChain.setTMScore(tmScore2);
			} catch (Exception e){
				e.printStackTrace();
				return afpChain;
			}

			//afpChain = CeMain.filterDuplicateAFPs(afpChain,calculator,ca1,ca2);
			//	prevAligs.add(afpChain);
			if ( displayJmol) {
				//if ( afpChain.getProbability() > 3.5) 
				showCurrentAlig(afpChain, ca1, ca2);
			}

			//System.out.println("We went through " + i + " iterations. Final TM score: " + afpChain.getTMScore());
			//afpChain.setAlgorithmName("CE-symmetry final result ");
			//StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			//jmol.evalString("draw l1 line 100 {0 0 0} (1:A.CA/1) ; draw l2 line 100 {0 0 0} (1:A.CA/2);" );
		} catch (Exception e){
			e.printStackTrace();
		}


		return afpChain;

	}

	public ConfigStrucAligParams getParameters() {
		return params;
	}

	public void setParameters(ConfigStrucAligParams parameters) {
		if ( ! (parameters instanceof CeParameters)){
			throw new IllegalArgumentException("Need to provide CeParameters, but provided " + parameters.getClass().getName());
		}

		params = (CeParameters) parameters;

	}

	public String getAlgorithmName() {
		return algorithmName;
	}

	public String getVersion() {
		return version;
	}


	/**
	 * Guesses the order of a symmetric alignment.
	 * 
	 * <p><strong>Details</strong><br/>
	 * Considers the distance (in number of residues) which a residue moves
	 * after undergoing <i>n</i> transforms by the alignment. If <i>n</i> corresponds
	 * to the intrinsic order of the alignment, this will be small. This algorithm
	 * tries increasing values of <i>n</i> and looks for abrupt decreases is the
	 * sum of squared distances. If none are found at <i>n</i><=8 (the maximum
	 * symmetry CE-Symm is likely to find), the alignment is reported as non-symmetric.
	 * @param afpChain A CE-symm alignment, where one protein is compared to itself
	 * @return The order of the alignment, or -1 if non-symmetric.
	 * @throws StructureException If afpChain is not one-to-one
	 */
	public static int getSymmetryOrder(AFPChain afpChain) throws StructureException {
		//maximum degree of rotational symmetry to consider
		final int maxSymmetry = 8;

		// Percentage change in RSSE required to improve score
		// Avoids reporting slight improvements in favor of lower order
		final float minimumMetricChange = 0.40f;

		Map<Integer,Integer> alignment = AlignmentTools.alignmentAsMap(afpChain);

		return AlignmentTools.getSymmetryOrder(alignment,
				new AlignmentTools.IdentityMap<Integer>(),
				maxSymmetry, minimumMetricChange);
	}


}
