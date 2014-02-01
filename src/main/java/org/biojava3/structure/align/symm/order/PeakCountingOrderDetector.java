package org.biojava3.structure.align.symm.order;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.util.Pair;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

/**
 * Determines order by smoothing and counting the number of peaks.
 * @author dmyersturnbull
 */
public class PeakCountingOrderDetector implements OrderDetector {

	private int maxOrder = 9;
	private double degreeSampling = 0.01;
	private double epsilon = 0.000001;
	private double bandwidth = 0.1;
	private int robustnessIterations = LoessInterpolator.DEFAULT_ROBUSTNESS_ITERS;
	private double loessAccuracy = LoessInterpolator.DEFAULT_ACCURACY;

	public PeakCountingOrderDetector(int maxOrder) {
		super();
		this.maxOrder = maxOrder;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);
			Pair<double[],double[]> pair = RotationOrderDetector.sampleRotations(ca, axis, degreeSampling);

			LoessInterpolator loess = new LoessInterpolator(bandwidth, robustnessIterations, loessAccuracy);

			double[] smoothed = loess.smooth(pair.getKey(), pair.getValue());

			for (double d : pair.getKey()) {
				System.out.print(StatUtils.formatD(d) + "\t");
			}
			System.out.println();

			for (double d : smoothed) {
				System.out.print(StatUtils.formatD(d) + "\t");
			}
			System.out.println();
			
			int nPeaks = countPeaks(smoothed, epsilon * Math.PI/180);
			
			/*
			 *  TODO Currently this isn't likely to handle order=1 well,
			 *  since C1 cases can easily have, say, exactly 5 peaks.
			 *  We will need to combine this with some smarter fitting method.
			 */

			return nPeaks; // for now
			
//			return nPeaks>maxOrder? 1 : nPeaks;

		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}

	}

	private int countPeaks(double[] values, double epsilon) {

		// TODO There's an off-by-1 error in odd cases
		// TODO I don't think we can actually use epsilon
		
		int nPeaks = 0;
		boolean previouslyIncreased = false;
		double previousValue = values[0]; // we can't have this count as previouslyIncreased
		for (double value : values) {
			if (previouslyIncreased && value < previousValue) {
				nPeaks++;
			}
			if (value > previousValue) {
				previouslyIncreased = true;
			} else {
				previouslyIncreased = false;
			}
			previousValue = value;
		}
		
		return nPeaks;

	}
	
	public void setMaxOrder(int maxOrder) {
		this.maxOrder = maxOrder;
	}

	public void setDegreeSampling(double degreeSampling) {
		this.degreeSampling = degreeSampling;
	}

	/**
	 * <em>After smoothing</em>, a threshold to count a peak.
	 * This should be very small but nonzero.
	 * This helps prevent spurious peaks from rounding errors from being counted.
	 * @param epsilon
	 */
	public void setEpsilon(double epsilon) {
		this.epsilon = epsilon;
	}

	public void setBandwidth(double bandwidth) {
		this.bandwidth = bandwidth;
	}

	public void setRobustnessIterations(int robustnessIterations) {
		this.robustnessIterations = robustnessIterations;
	}

	public void setLoessAccuracy(double loessAccuracy) {
		this.loessAccuracy = loessAccuracy;
	}

	public static void main(String[] args) throws Exception {

//		String name = "d1ijqa1"; // 6
		String name = "1TIM.A"; // 8
//		String name = "d1h70a_"; // 5

		// Perform alignment to determine axis
		Atom[] ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
		Atom[] ca2 = StructureTools.cloneCAArray(ca1);
		CeSymm ce = new CeSymm();
		AFPChain alignment = ce.align(ca1, ca2);

		// Search for orders up to 9
		PeakCountingOrderDetector detector = new PeakCountingOrderDetector(9);

		// Calculate order
		int order = detector.calculateOrder(alignment, ca1);
		System.out.println("Order: " + order);

	}
}
