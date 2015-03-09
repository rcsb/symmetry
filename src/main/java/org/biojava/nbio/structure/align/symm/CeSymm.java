package org.biojava.nbio.structure.align.symm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.AbstractStructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.ce.MatrixListener;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava.nbio.structure.align.symm.order.OrderDetector;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.align.symm.subunit.SubunitTools;
import org.biojava.nbio.structure.utils.SymmetryTools;

/**
 * Try to identify all possible symmetries by iterating recursively over all
 * results and disabling the diagonal of each previous result.
 * 
 * @author andreas
 * 
 * Modified Aleix Lafita: 09.03.2015
 * 
 */
public class CeSymm extends AbstractStructureAlignment implements
		MatrixListener, StructureAlignment {
	static final boolean debug = false;

	public static final String algorithmName = "jCE-symmetry";

	public static final String version = "1.0";

	private OrderDetector orderDetector = new SequenceFunctionOrderDetector(8,
			0.4f); // TODO finish

	AFPChain afpChain;

	Atom[] ca1;
	Atom[] ca2;
	int rows;
	int cols;
	CECalculator calculator;
	CESymmParameters params;
	// int loopCount ;
	
	public CeSymm() {
		super();
		params = new CESymmParameters();
	}

	public static void main(String[] args) {
		// Responsible for creating a CeMain instance
		CeSymmUserArgumentProcessor processor = new CeSymmUserArgumentProcessor(); 
		processor.process(args);
	}



	public static String toDBSearchResult(AFPChain afpChain) {
		StringBuffer str = new StringBuffer();

		str.append(afpChain.getName1());
		str.append("\t");
		str.append(afpChain.getName2());
		str.append("\t");
		str.append(String.format("%.2f", afpChain.getAlignScore()));
		str.append("\t");
		str.append(String.format("%.2f", afpChain.getProbability()));
		str.append("\t");
		str.append(String.format("%.2f", afpChain.getTotalRmsdOpt()));
		str.append("\t");
		str.append(afpChain.getCa1Length());
		str.append("\t");
		str.append(afpChain.getCa2Length());
		str.append("\t");
		str.append(afpChain.getCoverage1());
		str.append("\t");
		str.append(afpChain.getCoverage2());
		str.append("\t");
		str.append(String.format("%.2f", afpChain.getTMScore()));
		str.append("\t");
		str.append(afpChain.getOptLength());
		return str.toString();
	}

	public AFPChain indentifyAllSymmetries(String name1, String name2,
			AtomCache cache, int fragmentLength) throws StructureException,
			IOException {

		params = new CESymmParameters();

		params.setWinSize(fragmentLength);

		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);

		return align(ca1, ca2, params);

	}

	private static Matrix align(AFPChain afpChain, Atom[] ca1, Atom[] ca2,
			CESymmParameters params, Matrix origM, CECalculator calculator,
			int counter) throws StructureException {

		int fragmentLength = params.getWinSize();

		Atom[] ca2clone = SymmetryTools.cloneAtoms(ca2);

		int rows = ca1.length;
		int cols = ca2.length;

		// Matrix that tracks similarity of a fragment of length fragmentLength
		// starting a position i,j.

		int blankWindowSize = fragmentLength;
		if (origM == null) {
			// Build alignment ca1 to ca2-ca2

			afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);

			origM = SymmetryTools.blankOutPreviousAlignment(afpChain, ca2,
					rows, cols, calculator, null, blankWindowSize);

		} else {
			// we are doing an iteration on a previous alignment
			// mask the previous alignment
			origM = SymmetryTools.blankOutPreviousAlignment(afpChain, ca2,
					rows, cols, calculator, origM, blankWindowSize);
		}
		//System.out.println("origM: blankout correct...");

		Matrix clone = (Matrix) origM.clone();

		// that's the matrix to run the alignment on..
		calculator.setMatMatrix(clone.getArray());
		//System.out.println("origM: matrix set to calculator correct...");

		calculator.traceFragmentMatrix(afpChain, ca1, ca2clone);
		//System.out.println("origM: trace fragment matrix correct...");

		calculator.nextStep(afpChain, ca1, ca2clone);
		//System.out.println("origM: next step correct...");

		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);

		afpChain.setDistanceMatrix(origM);
		//System.out.println("origM: distance matrix correct...");
		
		return origM;

	}

	@Override
	public double[][] matrixInOptimizer(double[][] max) {

		return CECalculator.updateMatrixWithSequenceConservation(max, ca1, ca2,
				params);
	}

	@Override
	public boolean[][] initializeBreakFlag(boolean[][] breakFlag) {
		int fragmentLength = params.getWinSize();
		try {
			if (afpChain != null) {
				breakFlag = SymmetryTools.blankOutBreakFlag(afpChain, ca2,
						rows, cols, calculator, breakFlag, fragmentLength);
			}
		} catch (Exception e) {
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

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {

		if (params == null)
			params = new CESymmParameters();

		return align(ca1, ca2, params);
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2O, Object param)
			throws StructureException {
		if (!(param instanceof CESymmParameters))
			throw new IllegalArgumentException(
					"CE algorithm needs an object of call CESymmParameters as argument.");

		this.params = (CESymmParameters) param;

		this.ca1 = ca1;
		this.ca2 = ca2O;

		ca2 = StructureTools.duplicateCA2(ca2O);
		rows = ca1.length;
		cols = ca2.length;

		Matrix origM = null;

		AFPChain myAFP = new AFPChain();

		calculator = new CECalculator(params);
		calculator.addMatrixListener(this);

		int i = 0;

		while ((afpChain == null) && i < params.getMaxNrAlternatives()) {

			if (origM != null) {
				myAFP.setDistanceMatrix((Matrix) origM.clone());
			}
			origM = align(myAFP, ca1, ca2, params, origM, calculator, i);

			double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
			myAFP.setTMScore(tmScore2);

			i++;
		}
		afpChain = myAFP;

		try {
			afpChain = CeCPMain.postProcessAlignment(afpChain, ca1, ca2,
					calculator);
		} catch (Exception e) {
			e.printStackTrace();
			return afpChain;
		}

		if (params.isRefineResult()) {
			int order;
			try {
				order = orderDetector.calculateOrder(myAFP, ca1);
				afpChain = SymmRefiner.refineSymmetry(afpChain, ca1, ca2O,
						order);
			} catch (OrderDetectionFailedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		double tmScore2 = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
		afpChain.setTMScore(tmScore2);

		return afpChain;

	}
	
	/**
	 * New method that stores all the successive AFP alignment in the allAlignments list before each blackout.
	 * Guesses the order of symmetry by detecting a drop in the alignment length of 10%.
	 * 
	 * @author Aleix Lafita
	 * 
	 * @param ca1
	 * @param ca2O
	 * @param param
	 * @param allAlignments
	 * @throws StructureException
	 */
	public void align(Atom[] ca1, Atom[] ca2O, Object param, ArrayList<AFPChain> allAlignments)
			throws StructureException {
		if (!(param instanceof CESymmParameters))
			throw new IllegalArgumentException(
					"CE algorithm needs an object of call CESymmParameters as argument.");

		this.params = (CESymmParameters) param;

		this.ca1 = ca1;
		this.ca2 = ca2O;

		ca2 = StructureTools.duplicateCA2(ca2O);
		rows = ca1.length;
		cols = ca2.length;
		
		calculator = new CECalculator(params);
		calculator.addMatrixListener(this);

		Matrix origM = null;
		
		AFPChain myAFP = new AFPChain();

		Integer OptAlgnLenth = null;
		
		int i = 1;
		System.out.println("Start of the loop CeSym align...");

		while ((afpChain == null) && i < params.getMaxNrAlternatives()) {

			if (origM != null) {
				myAFP.setDistanceMatrix((Matrix) origM.clone());
			}
			origM = align(myAFP, ca1, ca2, params, origM, calculator, i+5);
			
			double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
			myAFP.setTMScore(tmScore2);
		
			//Clone the AFPChain
			AFPChain newAFP = (AFPChain) myAFP.clone();
			double tmScore3 = AFPChainScorer.getTMScore(newAFP, ca1, ca2);
			newAFP.setTMScore(tmScore3);
			
			//Post process the alignment
			try {
				System.out.println("Post process alignment "+i+"...");
				newAFP = CeCPMain.postProcessAlignment(newAFP, ca1, ca2,
						calculator);
			} catch (Exception e) {
				e.printStackTrace();
				allAlignments.add(newAFP);
			}
			
			//If it is the first alignment set the optimal length
			if (OptAlgnLenth==null){
				OptAlgnLenth = newAFP.getOptLength();
			}
			//If not check for a drop in the alignment length and break the loop
			else if (newAFP.getOptLength() < (OptAlgnLenth-OptAlgnLenth/10)){
				System.out.println("Order of symmetry detected: "+i);
				System.out.println("Optimal Alignment Lenth: "+OptAlgnLenth+", Last alignment length: "+newAFP.getOptLength());
				break;
			}
			//Add the alignment to the allAlignments list otherwise
			allAlignments.add(newAFP);
			System.out.println("Alignment "+i+" completed...");
			
			i++;
		}
		
		System.out.println("CeSym align completed...");
	}

	/**
	 * New method that creates a multiple block AFP alignment corresponding to the subunits of symmetry.
	 * Guesses the order of symmetry by detecting a drop in the alignment length of 10% (and score to implement).
	 * 
	 * @author Aleix Lafita
	 * 
	 */
	public AFPChain alignMultiple(Atom[] ca1, Atom[] ca2O, Object param)
			throws StructureException {
		if (!(param instanceof CESymmParameters))
			throw new IllegalArgumentException(
					"CE algorithm needs an object of call CESymmParameters as argument.");

		this.params = (CESymmParameters) param;

		this.ca1 = ca1;
		this.ca2 = ca2O;

		ca2 = StructureTools.duplicateCA2(ca2O);
		rows = ca1.length;
		cols = ca2.length;
		
		calculator = new CECalculator(params);
		calculator.addMatrixListener(this);

		Matrix origM = null;
		
		AFPChain myAFP = new AFPChain();
		List<int[][][]> allAlignments = new ArrayList<int[][][]>();

		Integer OptAlgnLen = null;
		
		int i = 1;
		System.out.println("Start of the loop CeSym align...");

		while (afpChain == null) {

			if (origM != null) {
				myAFP.setDistanceMatrix((Matrix) origM.clone());
			}
			System.out.println("Align matrix number "+i+"...");
			origM = align(myAFP, ca1, ca2, params, origM, calculator, i+5);
			
			System.out.println("Get TM score number "+i+"...");
			double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
			myAFP.setTMScore(tmScore2);
			
			System.out.println("myAFP alignment "+i+" has length "+myAFP.getOptAln()[0][0].length+"...");
		
			//Clone the AFPChain
			AFPChain newAFP = (AFPChain) myAFP.clone();
			System.out.println("AFPChain cloned number "+i+"...");
			
			//Post process the alignment
			try {
				System.out.println("Post process alignment "+i+"...");
				newAFP = CeCPMain.postProcessAlignment(newAFP, ca1, ca2,
						calculator);
			} catch (Exception e) {
				System.out.println("Post process alignment in CeSym align failed...");
				e.printStackTrace();
				allAlignments.add(newAFP.getOptAln().clone());
			}
			
			System.out.println("newAFP alignment "+i+" has length "+newAFP.getOptAln()[0][0].length+"...");
			
			//NEEDED...?
			System.out.println("Set TM score number "+i+"...");
			double tmScore3 = AFPChainScorer.getTMScore(newAFP, ca1, ca2);
			newAFP.setTMScore(tmScore3);
			
			//If it is the first alignment set the optimal length
			if (OptAlgnLen==null){
				OptAlgnLen = newAFP.getOptLength();
			}
			//If not check for a drop in the alignment length and break the loop
			else if (newAFP.getOptLength() < (OptAlgnLen-OptAlgnLen/10)){
				System.out.println("Order of symmetry detected: "+i);
				System.out.println("Optimal Alignment Lenth: "+OptAlgnLen+", Last alignment length: "+newAFP.getOptLength());
				break;
			}
			//Add the alignment to the allAlignments list otherwise
			allAlignments.add(newAFP.getOptAln().clone());
			System.out.println("Alignment "+i+" completed...");
			
			i++;
			
		}
		
		System.out.println("CeSym align completed...");
		
		afpChain = SubunitTools.refinedAFP(allAlignments, ca1);
		return afpChain;
		
	}
	
	public AFPChain alignMultiple(Atom[] ca1, Atom[] ca2) throws StructureException {

		if (params == null)
			params = new CESymmParameters();

		return alignMultiple(ca1, ca2, params);
	}

	@Override
	public ConfigStrucAligParams getParameters() {
		return params;
	}

	@Override
	public void setParameters(ConfigStrucAligParams parameters) {
		if (!(parameters instanceof CESymmParameters)) {
			throw new IllegalArgumentException(
					"Need to provide CESymmParameters, but provided "
							+ parameters.getClass().getName());
		}

		params = (CESymmParameters) parameters;
	}

	@Override
	public String getAlgorithmName() {
		return algorithmName;
	}

	@Override
	public String getVersion() {
		return version;
	}

	public OrderDetector getOrderDetector() {
		return orderDetector;
	}

	public void setOrderDetector(OrderDetector orderDetector) {
		this.orderDetector = orderDetector;
	}

	public static boolean isSignificant(AFPChain afpChain,OrderDetector orderDetector, Atom[] ca1) throws StructureException {

		// TM-score cutoff
		if (afpChain.getTMScore() < 0.4) return false;

		// sequence-function order cutoff
		int order = 1;
			try {
				order = orderDetector.calculateOrder(afpChain, ca1);
			} catch (OrderDetectionFailedException e) {
				e.printStackTrace();
				// try the other method
			}

		if (order > 1) return true;

		// angle order cutoff
		RotationAxis rot = new RotationAxis(afpChain);
		order = rot.guessOrderFromAngle(1.0 * Calc.radiansPerDegree, 8);
		if (order > 1) return true;
		// asymmetric
		return false;
	}
	
	public boolean isSignificant() throws StructureException {
		return CeSymm.isSignificant(this.afpChain,this.orderDetector,this.ca1);
	}
}
