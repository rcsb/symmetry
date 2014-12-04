package org.biojava3.structure.align.symm;

import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.AbstractStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.ce.MatrixListener;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava3.structure.align.symm.order.OrderDetector;
import org.biojava3.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava3.structure.utils.SymmetryTools;

/**
 * Try to identify all possible symmetries by iterating resursively over all
 * results and disabling the diagonal of each previous result.
 * 
 * @author andreas
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

		Matrix clone = (Matrix) origM.clone();

		// that's the matrix to run the alignment on..
		calculator.setMatMatrix(clone.getArray());

		calculator.traceFragmentMatrix(afpChain, ca1, ca2clone);

		calculator.nextStep(afpChain, ca1, ca2clone);

		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);

		afpChain.setDistanceMatrix(origM);
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

			afpChain = myAFP;
			if (origM != null) {
				myAFP.setDistanceMatrix((Matrix) origM.clone());
			}
			origM = align(myAFP, ca1, ca2, params, origM, calculator, i);

			double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
			myAFP.setTMScore(tmScore2);

			i++;
		}

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
