package org.biojava3.structure.quaternary.core;

import java.util.Arrays;

public class QuatSymmetryParameters {
	private int minimumSequenceLength = 20;
	private double[] sequenceIdentityThresholds = {0.0, 0.95};
	private double sequencePseudoSymmetryThreshold = 0.95;
	private double alignmentFractionThreshold = 0.9;
	private double rmsdThreshold = 7.0;
	private double angleThreshold = 10.0;
	private int maximumLocalCombinations = 100000;
	private boolean localSymmetry = true;
	private boolean verbose = false;
	private static final String n = System.getProperty("line.separator");
	
	/**
	 * @return the minimumSequenceLength
	 */
	public int getMinimumSequenceLength() {
		return minimumSequenceLength;
	}
	/**
	 * @param minimumSequenceLength the minimumSequenceLength to set
	 */
	public void setMinimumSequenceLength(int minimumSequenceLength) {
		this.minimumSequenceLength = minimumSequenceLength;
	}
	/**
	 * @return the sequenceIdentityThreshold
	 */
	
	/**
	 * @return the sequenceIdentityThreshold
	 */
	public double[] getSequenceIdentityThresholds() {
		return sequenceIdentityThresholds;
	}
	
	/**
	 * @param sequenceIdentityThresholds the sequenceIdentityThresholds to set
	 */
	public void setSequenceIdentityThresholds(double[] sequenceIdentityThresholds) {
		this.sequenceIdentityThresholds = sequenceIdentityThresholds;
	}
	
	/**
	 * @return the alignmentFractionThreshold
	 */
	public double getAlignmentFractionThreshold() {
		return alignmentFractionThreshold;
	}
	/**
	 * @param alignmentFractionThreshold the alignmentFractionThreshold to set
	 */
	public void setAlignmentFractionThreshold(double alignmentFractionThreshold) {
		this.alignmentFractionThreshold = alignmentFractionThreshold;
	}
	/**
	 * @return the rmsdThreshold
	 */
	public double getRmsdThreshold() {
		return rmsdThreshold;
	}
	/**
	 * @param rmsdThreshold the rmsdThreshold to set
	 */
	public void setRmsdThreshold(double rmsdThreshold) {
		this.rmsdThreshold = rmsdThreshold;
	}
	public double getAngleThreshold() {
		return angleThreshold;
	}
	public void setAngleThreshold(double angleThreshold) {
		this.angleThreshold = angleThreshold;
	}
	
	public double getSequencePseudoSymmetryThreshold() {
		return sequencePseudoSymmetryThreshold;
	}
	
	public void setSequencePseudoSymmetryThreshold(
			double sequencePseudoSymmetryThreshold) {
		this.sequencePseudoSymmetryThreshold = sequencePseudoSymmetryThreshold;
	}
	
	public int getMaximumLocalCombinations() {
		return maximumLocalCombinations;
	}
	public void setMaximumLocalCombinations(int maximumLocalCombinations) {
		this.maximumLocalCombinations = maximumLocalCombinations;
	}
	public boolean isLocalSymmetry() {
		return localSymmetry;
	}
	public void setLocalSymmetry(boolean localSymmetry) {
		this.localSymmetry = localSymmetry;
	}
	public boolean isVerbose() {
		return verbose;
	}
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}
	
	public String toString() {
		StringBuilder s = new StringBuilder();
		s.append("Minimum protein sequence length   : ");
		s.append(minimumSequenceLength);
		s.append(n);
		s.append("Sequence identity thresholds      : ");
		s.append(Arrays.toString(sequenceIdentityThresholds));
		s.append(n);
		s.append("Sequence pseudosymmetry threshold : ");
		s.append(sequencePseudoSymmetryThreshold);
		s.append(n);
		s.append("Alignment fraction threshold      : ");
		s.append(alignmentFractionThreshold);
		s.append(n);
		s.append("Angle threshold                   : ");
		s.append(angleThreshold);
		s.append(n);	
		s.append("Symmetry RMSD threshold           : ");
		s.append(rmsdThreshold);
		s.append(n);
		s.append("Local symmetry                    : ");
		s.append(localSymmetry);
		s.append(n);
		s.append("Verbose                           : ");
		s.append(verbose);
		s.append(n);
		return s.toString();
	}
}
