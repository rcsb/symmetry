package org.biojava3.structure.quaternary.core;

public class QuatSymmetryParameters {
	private int minimumSequenceLength = 24;
	private boolean structuralAlignmentOnly = false;
	private double sequenceIdentityThreshold = 0.30;
	private double sequencePseudoSymmetryThreshold = 0.95;
	private double alignmentFractionThreshold = 0.9;
	private double rmsdThreshold = 5.0;
	
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
	
	public String toString() {
		StringBuilder s = new StringBuilder();
		s.append("Minimum protein sequence length  : " + minimumSequenceLength  +"/n");
		s.append("Sequence identity threshold      : " + sequenceIdentityThreshold + "/n");
		s.append("Sequence pseudosymmetry threshold: " + sequencePseudoSymmetryThreshold + "/n");
		s.append("Alignment fraction threshold     : " + alignmentFractionThreshold + "/n");
		s.append("Symmetry RMSD threshold          : " + rmsdThreshold + "/n");
		return s.toString();
	}
}
