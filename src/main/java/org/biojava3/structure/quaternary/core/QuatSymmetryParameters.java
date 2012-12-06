package org.biojava3.structure.quaternary.core;

public class QuatSymmetryParameters {
	private int minimumSequenceLength = 24;
	private boolean structuralAlignmentOnly = false;
	private double sequenceIdentityThreshold = 0.30;
	private double sequencePseudoSymmetryThreshold = 0.95;
	private double alignmentFractionThreshold = 0.9;
	private double rmsdThreshold = 5.0;
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
	public double getSequenceIdentityThreshold() {
		return sequenceIdentityThreshold;
	}
	/**
	 * @param sequenceIdentityThreshold the sequenceIdentityThreshold to set
	 */
	public void setSequenceIdentityThreshold(double sequenceIdentityThreshold) {
		this.sequenceIdentityThreshold = sequenceIdentityThreshold;
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
	/**
	 * @return the structuralAlignmentOnly
	 */
	public boolean isStructuralAlignmentOnly() {
		return structuralAlignmentOnly;
	}
	
// Reserved for the future. This feature is currently not supported
//	/**
//	 * @param structuralAlignmentOnly the structuralAlignmentOnly to set
//	 */
//	public void setStructuralAlignmentOnly(boolean structuralAlignmentOnly) {
//		this.structuralAlignmentOnly = structuralAlignmentOnly;
//	}
	
	public double getSequencePseudoSymmetryThreshold() {
		return sequencePseudoSymmetryThreshold;
	}
	
	public void setSequencePseudoSymmetryThreshold(
			double sequencePseudoSymmetryThreshold) {
		this.sequencePseudoSymmetryThreshold = sequencePseudoSymmetryThreshold;
	}
	
	public boolean isVerbose() {
		return verbose;
	}
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}
	
	public String toString() {
		StringBuilder s = new StringBuilder();
		s.append("Minimum protein sequence length  : ");
		s.append(minimumSequenceLength);
		s.append(n);
		s.append("Sequence identity threshold      : ");
		s.append(sequenceIdentityThreshold);
		s.append(n);
		s.append("Sequence pseudosymmetry threshold: ");
		s.append(sequencePseudoSymmetryThreshold);
		s.append(n);
		s.append("Alignment fraction threshold     : ");
		s.append(alignmentFractionThreshold);
		s.append(n);
		s.append("Symmetry RMSD threshold          : ");
		s.append(rmsdThreshold);
		s.append(n);
		s.append("Verbose                          : ");
		s.append(verbose);
		s.append(n);
		return s.toString();
	}
}
