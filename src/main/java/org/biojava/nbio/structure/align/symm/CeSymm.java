package org.biojava.nbio.structure.align.symm;

import java.awt.Color;
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
import org.biojava.nbio.structure.align.symm.CESymmParameters.SubunitColors;
import org.biojava.nbio.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava.nbio.structure.align.symm.order.OrderDetector;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.align.symm.refine.MultipleAlignRefiner;
import org.biojava.nbio.structure.align.symm.refine.Refiner;
import org.biojava.nbio.structure.align.symm.refine.RefinerFailedException;
import org.biojava.nbio.structure.align.symm.refine.SingleAlignRefiner;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.utils.SymmetryTools;
import org.jcolorbrewer.ColorBrewer;

/**
 * Try to identify all possible symmetries by iterating recursively over all
 * results and disabling the diagonal of each previous result.
 * 
 * @author andreas
 * 
 * 
 */
public class CeSymm extends AbstractStructureAlignment implements
		MatrixListener, StructureAlignment {
	static final boolean debug = false;

	public static final String algorithmName = "jCE-symmetry";

	public static final String version = "1.0";
	
	//The order and refinement options are controlled by CESymmParameters
	private OrderDetector orderDetector = new SequenceFunctionOrderDetector(8, 0.4f);
	private Refiner refiner = null;

	AFPChain afpChain;
	AFPChain[] afpAlignments;

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
	public AFPChain align(Atom[] ca10, Atom[] ca2O, Object param)
			throws StructureException {
		if (!(param instanceof CESymmParameters))
			throw new IllegalArgumentException(
					"CE algorithm needs an object of call CESymmParameters as argument.");

		this.params = (CESymmParameters) param;

		ca1 = ca10;
		ca2 = StructureTools.duplicateCA2(ca2O);
		rows = ca1.length;
		cols = ca2.length;

		Matrix origM = null;

		AFPChain myAFP = new AFPChain();
		List<AFPChain> allAlignments = new ArrayList<AFPChain>();
		
		//The two variables that keep track of the goodness of the alignment compared to the optimal
		int optAlgnLen = 0;
		double optTMscore = 0;

		calculator = new CECalculator(params);
		calculator.addMatrixListener(this);
		
		//Set multiple to true if multiple alignments are needed
		boolean multiple = (params.getRefineMethod() == RefineMethod.MULTIPLE);

		int i = 0;

		do {
			//System.out.print("Alignment number: "+i);

			if (origM != null) {
				myAFP.setDistanceMatrix((Matrix) origM.clone());
			}
			origM = align(myAFP, ca1, ca2, params, origM, calculator, i);

			double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
			myAFP.setTMScore(tmScore2);
			
			//Clone the AFPChain
			AFPChain newAFP = (AFPChain) myAFP.clone();
			
			//Post process the alignment
			try {
				newAFP = CeCPMain.postProcessAlignment(newAFP, ca1, ca2,
						calculator);
			} catch (Exception e) {
				e.printStackTrace();
				allAlignments.add(newAFP);
			}
			
			//Calculate and set the TM score for the newAFP alignment before adding it to the list
			double tmScore3 = AFPChainScorer.getTMScore(newAFP, ca1, ca2);
			newAFP.setTMScore(tmScore3);
			
			//Print alignment length and TM score of the current alignment
			//System.out.println("Alignment "+(i+1)+" length: "+newAFP.getOptLength());
			//System.out.println("Alignment "+(i+1)+" score: "+newAFP.getTMScore());
			
			//If it is the first alignment set the optimal length
			if (i==0){
				optAlgnLen = newAFP.getOptLength();
				optTMscore = newAFP.getTMScore();
				//Threshold for symmetry detection TBD
				if (optTMscore < 0.2){
					System.out.println("No symmetry detected");
					return newAFP;
				}
			}
			//If not check for a drop (of 50%) in the alignment length or score and break the loop
			else if ((newAFP.getOptLength() < (optAlgnLen-optAlgnLen/2) || newAFP.getTMScore() < (optTMscore-optTMscore/2))){
				//System.out.println("Optimal alignment length: "+optAlgnLen+", Last alignment length: "+newAFP.getOptLength());
				//System.out.println("Optimal alignment TM score: "+optTMscore+", Last alignment TM score: "+newAFP.getTMScore());
				break;
			}
			//Add the alignment to the allAlignments list
			allAlignments.add(newAFP);
			//System.out.println("Alignment "+(i+1)+" completed...");
			
			i++;
		} while (i < params.getMaxNrSubunits() && multiple);
		
		int order = params.getMaxNrSubunits();
		//Save the results to the AFPChain variables
		afpChain = allAlignments.get(0);
		afpAlignments = new AFPChain[allAlignments.size()];
		for (int k=0; k<allAlignments.size(); k++){
			afpAlignments[k] = allAlignments.get(k);
		}
		
		//Refinement options
		if (params.getRefineMethod() == RefineMethod.MULTIPLE){
			order = afpAlignments.length+1;
			System.out.println("Order of symmetry: "+(order));
			refiner = new MultipleAlignRefiner();
			try {
				afpChain = refiner.refine(afpAlignments, ca1, ca2, order);
			} catch (RefinerFailedException e) {
				e.printStackTrace();
			}
		}
		else if (params.getRefineMethod() == RefineMethod.SINGLE){
			refiner = new SingleAlignRefiner();
			//Calculate order
			try {
				order = orderDetector.calculateOrder(afpChain, ca1);
				System.out.println("Order of symmetry: "+(order));
			} catch (OrderDetectionFailedException e) {
				e.printStackTrace();
			}
			//Refine the AFPChain
			try {
				afpChain = refiner.refine(afpAlignments, ca1, ca2, order);
			} catch (RefinerFailedException e1) {
				e1.printStackTrace();
			}
		}
		
		//Coloring options
		if (params.getSubunitColors() == SubunitColors.COLOR_SET){
			Color[] colors = ColorBrewer.Set1.getColorPalette(afpChain.getBlockNum());
			afpChain.setBlockColors(colors);
		}
		else if (params.getSubunitColors() == SubunitColors.SPECTRAL){
			Color[] colors = ColorBrewer.Spectral.getColorPalette(afpChain.getBlockNum());
			afpChain.setBlockColors(colors);
		}

		double tmScore2 = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
		afpChain.setTMScore(tmScore2);
		
		//System.out.println("CeSymm alignment completed...");

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
	
	public Refiner getRefiner() {
		return refiner;
	}
	
	public void setRefiner(Refiner refiner) {
		this.refiner = refiner;
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
