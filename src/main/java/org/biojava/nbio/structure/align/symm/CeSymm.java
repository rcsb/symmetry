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
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.order.AngleOrderDetector;
import org.biojava.nbio.structure.align.symm.order.MultiMethodOrderDetector;
import org.biojava.nbio.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava.nbio.structure.align.symm.order.OrderDetector;
import org.biojava.nbio.structure.align.symm.order.PeakCountingOrderDetector;
import org.biojava.nbio.structure.align.symm.order.RotationOrderDetector;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.align.symm.refine.SymmOptimizer;
import org.biojava.nbio.structure.align.symm.refine.MultipleRefiner;
import org.biojava.nbio.structure.align.symm.refine.Refiner;
import org.biojava.nbio.structure.align.symm.refine.RefinerFailedException;
import org.biojava.nbio.structure.align.symm.refine.SingleRefiner;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Try to identify the possible symmetries and repeats in a structure by running an alignment of the structure
 * against itself disabling the diagonal of the identity alignment. Iterating recursively over all results
 * and disabling the diagonal of each previous result can also be done with the current implementation.
 * 
 * @author andreas
 * 
 */
public class CeSymm extends AbstractStructureAlignment implements MatrixListener, StructureAlignment {

	private static final boolean debug = false;
	public static final String algorithmName = "jCE-symmetry";
	public static final String version = "1.0";
	private static final Logger logger = LoggerFactory.getLogger(CeSymm.class);
	private static final double symmetryThreshold = 0.4;

	private AFPChain afpChain;
	private List<AFPChain> afpAlignments;

	private Atom[] ca1;
	private Atom[] ca2;
	private int rows;
	private int cols;
	
	private CECalculator calculator;
	private CESymmParameters params;
	
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

	private static Matrix align(AFPChain afpChain, Atom[] ca1, Atom[] ca2, CESymmParameters params, 
			Matrix origM, CECalculator calculator, int counter) throws StructureException {

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
		
		final Matrix origMfinal = (Matrix) origM.clone();
		//Add a matrix listener to keep the blacked zones in max.
		calculator.addMatrixListener(new MatrixListener() {

			@Override
			public double[][] matrixInOptimizer(double[][] max) {
				
				//Check every entry of origM for blacked out regions
				for (int i=0; i<max.length; i++){
					for (int j=0; j<max[i].length; j++){
						if (origMfinal.getArray()[i][j]>1e9){
							max[i][j] = -origMfinal.getArray()[i][j];
						}
					}
				}
				return max;
			}

			@Override
			public boolean[][] initializeBreakFlag(boolean[][] brkFlag) {
				
				return brkFlag;
			}
		});
		
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
	public AFPChain align(Atom[] ca10, Atom[] ca2O, Object param) throws StructureException {
		
	//STEP 0: prepare all the information for the symmetry alignment
		if (!(param instanceof CESymmParameters))
			throw new IllegalArgumentException("CE algorithm needs an object of call CESymmParameters as argument.");

		this.params = (CESymmParameters) param;

		ca1 = ca10;
		ca2 = StructureTools.duplicateCA2(ca2O);
		rows = ca1.length;
		cols = ca2.length;
		
		if( rows == 0 || cols == 0) {
			logger.warn("Aligning empty structure");
			throw new StructureException("Aligning empty structure");
		}

		Matrix origM = null;
		AFPChain myAFP = new AFPChain();
		afpAlignments = new ArrayList<AFPChain>();
		
		calculator = new CECalculator(params);
		calculator.addMatrixListener(this);
		
		//Set multiple to true if multiple alignments are needed for refinement
		boolean multiple = (params.getRefineMethod() == RefineMethod.MULTIPLE || params.getRefineMethod() == RefineMethod.MULTIPLE_OPTIMIZE);
	
	//STEP 1: perform the raw symmetry alignment
		int i = 0;
		do {
			
			if (origM != null) {
				myAFP.setDistanceMatrix((Matrix) origM.clone());
			}
			origM = align(myAFP, ca1, ca2, params, origM, calculator, i);

			double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
			myAFP.setTMScore(tmScore2);
			
			//Clone the AFPChain
			AFPChain newAFP = (AFPChain) myAFP.clone();
			
			//Post process the alignment
			newAFP = CeCPMain.postProcessAlignment(newAFP, ca1, ca2, calculator);
			
			//Calculate and set the TM score for the newAFP alignment
			double tmScore3 = AFPChainScorer.getTMScore(newAFP, ca1, ca2);
			newAFP.setTMScore(tmScore3);
			
			if (debug) logger.info("Alignment "+(i+1)+" score: "+newAFP.getTMScore());
			
			//Determine if the alignment is significant to do more alignment iterations
			if (!isSignificant(newAFP,ca1)){
				if(debug) logger.info("Not symmetric alignment with TM score: "+newAFP.getTMScore());
				//If it is the first alignment save it anyway and try to optimize it
				if (i==0) afpAlignments.add(newAFP);
				break;
			}
			//If it is a symmetric alignment add it to the allAlignments list
			afpAlignments.add(newAFP);
			
			i++;
			
		} while (i < params.getMaxSymmOrder() && multiple);
		
		//Save the results to the CeSymm member variables
		afpChain = afpAlignments.get(0);
		
	//STEP 2: calculate the order of symmetry / number of internal repeats
		int order = 1;
		OrderDetector orderDetector = null;
		switch (params.getOrderDetectorMethod()) {
		case SEQUENCE_FUNCTION: 
			orderDetector = new SequenceFunctionOrderDetector(params.getMaxSymmOrder(), 0.4f);
			break;
		case MULTI_METHOD:
			orderDetector = new MultiMethodOrderDetector(100, 1.0); //TODO parameters?
			break;
		case ANGLE:
			orderDetector = new AngleOrderDetector(params.getMaxSymmOrder(), 1.0);
			break;
		case PEAK_COUNTING:
			orderDetector = new PeakCountingOrderDetector(params.getMaxSymmOrder());
			break;
		case ROTATION:
			orderDetector = new RotationOrderDetector(params.getMaxSymmOrder());
			break;
		}
		try {
			order = orderDetector.calculateOrder(afpChain, ca1);
		} catch (OrderDetectionFailedException e) {
			e.printStackTrace();
		}
		
	//STEP 3: symmetry refinement, apply consistency in the subunit residues
		Refiner refiner = null;
		switch (params.getRefineMethod()){
		case MULTIPLE:
			refiner = new MultipleRefiner();
			break;
		case MULTIPLE_OPTIMIZE:
			refiner = new MultipleRefiner();
			break;
		case SINGLE:
			refiner = new SingleRefiner();
			break;
		case SINGLE_OPTIMIZE:
			refiner = new SingleRefiner();
			break;
		case NOT_REFINED:
			return afpChain;
		}
		try {
			afpChain = refiner.refine(afpAlignments, ca1, ca2, order);
		} catch (RefinerFailedException e) {
			e.printStackTrace();
		}
		
	//STEP 4: symmetry alignment optimization
		if (params.getRefineMethod() == RefineMethod.MULTIPLE_OPTIMIZE || params.getRefineMethod() == RefineMethod.SINGLE_OPTIMIZE){
			SymmOptimizer optimizer = new SymmOptimizer();
			try {
				afpChain = optimizer.optimize(afpChain, ca1, ca2, order);
			} catch (RefinerFailedException e) {
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
	
	public static boolean isSignificant(AFPChain afpChain, Atom[] ca1) throws StructureException {

		// TM-score cutoff
		if (afpChain.getTMScore() < symmetryThreshold) return false;

		// sequence-function order cutoff
		int order = 1;
			try {
				OrderDetector orderDetector = new SequenceFunctionOrderDetector(8, 0.4f);
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
		return CeSymm.isSignificant(this.afpChain,this.ca1);
	}
}
