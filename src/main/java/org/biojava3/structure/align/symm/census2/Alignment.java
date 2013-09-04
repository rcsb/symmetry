/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2013-02-22
 *
 */
package org.biojava3.structure.align.symm.census2;

import java.io.Serializable;

import org.biojava.bio.structure.align.model.AFPChain;

/**
 * The results of a symmetry alignment. Includes information for CE-Symm or SymD.
 * 
 * @author dmyerstu
 */
public class Alignment implements Serializable {

	private static final long serialVersionUID = 5560951230435020821L;

	private Integer alignLength;

	private Float alignScore;

	private Float tScore;

	private Integer block1Length;

	private Integer block2Length;

	private Integer coverage;

	private Integer gapLength;

	private Float identity;

	private Integer initialShift;

	private Integer nNonSelfAligned;

	private Float rmsd;

	private Float similarity;

	private Float symDZScore;
	
	private Float symDTMScore;

	private Float tmScore;

	private Float zScore;

	public Alignment() {

	}

	public Alignment(AFPChain afpChain) {
		if (afpChain != null) {
			identity = (float) afpChain.getIdentity();
			similarity = (float) afpChain.getSimilarity();
			coverage = afpChain.getCoverage1();
			zScore = (float) afpChain.getProbability();
			rmsd = (float) afpChain.getTotalRmsdOpt();
			tmScore = (float) afpChain.getTMScore();
			alignLength = afpChain.getAlnLength();
			gapLength = afpChain.getGapLen();
			alignScore = (float) afpChain.getAlignScore();
			block1Length = afpChain.getBlockResSize()[0];
			if (afpChain.getBlockResSize().length > 1) {
				block2Length = afpChain.getBlockResSize()[1];
			} else {
				block2Length = 0;
			}
		}
	}

	public Alignment(Integer alignLength, Integer coverage, Float alignScore, Integer gapLength, Float identity,
			Float rmsd, Float similarity, Float tmScore, Float zScore, Integer block1Length, Integer block2Length) {
		this.alignLength = alignLength;
		this.coverage = coverage;
		this.alignScore = alignScore;
		this.gapLength = gapLength;
		this.identity = identity;
		this.rmsd = rmsd;
		this.similarity = similarity;
		this.tmScore = tmScore;
		this.zScore = zScore;
		this.block1Length = block1Length;
		this.block2Length = block2Length;
	}

	public Integer getAlignLength() {
		return alignLength;
	}

	public Float getAlignScore() {
		return alignScore;
	}

	public Float getAlternateTm() {
		return tScore;
	}

	public Integer getBlock1Length() {
		return block1Length;
	}

	public Integer getBlock2Length() {
		return block2Length;
	}

	public Integer getCoverage() {
		return coverage;
	}

	public Integer getGapLength() {
		return gapLength;
	}

	public Float getIdentity() {
		return identity;
	}

	public Integer getInitialShift() {
		return initialShift;
	}

	public Integer getnNonSelfAligned() {
		return nNonSelfAligned;
	}

	public Float getRmsd() {
		return rmsd;
	}

	public Float getSimilarity() {
		return similarity;
	}

	public Float getTmScore() {
		return tmScore;
	}

	public Float getzScore() {
		return zScore;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (alignLength == null ? 0 : alignLength.hashCode());
		result = prime * result + (alignScore == null ? 0 : alignScore.hashCode());
		result = prime * result + (tScore == null ? 0 : tScore.hashCode());
		result = prime * result + (block1Length == null ? 0 : block1Length.hashCode());
		result = prime * result + (block2Length == null ? 0 : block2Length.hashCode());
		result = prime * result + (coverage == null ? 0 : coverage.hashCode());
		result = prime * result + (gapLength == null ? 0 : gapLength.hashCode());
		result = prime * result + (identity == null ? 0 : identity.hashCode());
		result = prime * result + (initialShift == null ? 0 : initialShift.hashCode());
		result = prime * result + (nNonSelfAligned == null ? 0 : nNonSelfAligned.hashCode());
		result = prime * result + (rmsd == null ? 0 : rmsd.hashCode());
		result = prime * result + (similarity == null ? 0 : similarity.hashCode());
		result = prime * result + (tmScore == null ? 0 : tmScore.hashCode());
		result = prime * result + (symDTMScore == null ? 0 : symDTMScore.hashCode());
		result = prime * result + (zScore == null ? 0 : zScore.hashCode());
		return result;
	}

	public void setAlignLength(Integer alignLength) {
		this.alignLength = alignLength;
	}

	public void setAlignScore(Float alignScore) {
		this.alignScore = alignScore;
	}

	public void setAlternateTm(Float alternateTm) {
		this.tScore = alternateTm;
	}

	public void setBlock1Length(Integer block1Length) {
		this.block1Length = block1Length;
	}

	public void setBlock2Length(Integer block2Length) {
		this.block2Length = block2Length;
	}

	public void setCoverage(Integer coverage) {
		this.coverage = coverage;
	}

	public void setGapLength(Integer gapLength) {
		this.gapLength = gapLength;
	}

	public void setIdentity(Float identity) {
		this.identity = identity;
	}

	public void setInitialShift(Integer initialShift) {
		this.initialShift = initialShift;
	}

	public void setnNonSelfAligned(Integer nNonSelfAligned) {
		this.nNonSelfAligned = nNonSelfAligned;
	}

	public void setRmsd(Float rmsd) {
		this.rmsd = rmsd;
	}

	public void setSimilarity(Float similarity) {
		this.similarity = similarity;
	}

	public void setTmScore(Float tmScore) {
		this.tmScore = tmScore;
	}

	public void setzScore(Float zScore) {
		this.zScore = zScore;
	}

	public Float gettScore() {
		return tScore;
	}

	public void settScore(Float tScore) {
		this.tScore = tScore;
	}

	public Float getSymDZScore() {
		return symDZScore;
	}

	public void setSymDZScore(Float symDZScore) {
		this.symDZScore = symDZScore;
	}

	public Float getSymDTMScore() {
		return symDTMScore;
	}

	public void setSymDTMScore(Float symDTMScore) {
		this.symDTMScore = symDTMScore;
	}

	@Override
	public String toString() {
		return "Alignment [alignLength=" + alignLength + ", alignScore=" + alignScore + ", block1Length="
				+ block1Length + ", block2Length=" + block2Length + ", coverage=" + coverage + ", gapLength="
				+ gapLength + ", identity=" + identity + ", rmsd=" + rmsd + ", similarity=" + similarity + ", tmScore="
				+ tmScore + ", alternateTm=" + tScore + ", zScore=" + zScore + ", initialShift=" + initialShift
				+ ", nNonSelfAligned=" + nNonSelfAligned + ", tmpr=" + symDTMScore + "]";
	}

}
