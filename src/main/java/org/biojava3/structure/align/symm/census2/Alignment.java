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

	private Float alternateTm;

	private Integer block1Length;

	private Integer block2Length;

	private Integer coverage;

	private Integer gapLength;

	private Float identity;

	private Integer initialShift;

	private Integer nNonSelfAligned;

	private Float rmsd;

	private Float similarity;

	private Float tmpr;

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

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Alignment other = (Alignment) obj;
		if (alignLength == null) {
			if (other.alignLength != null)
				return false;
		} else if (!alignLength.equals(other.alignLength))
			return false;
		if (alignScore == null) {
			if (other.alignScore != null)
				return false;
		} else if (!alignScore.equals(other.alignScore))
			return false;
		if (alternateTm == null) {
			if (other.alternateTm != null)
				return false;
		} else if (!alternateTm.equals(other.alternateTm))
			return false;
		if (block1Length == null) {
			if (other.block1Length != null)
				return false;
		} else if (!block1Length.equals(other.block1Length))
			return false;
		if (block2Length == null) {
			if (other.block2Length != null)
				return false;
		} else if (!block2Length.equals(other.block2Length))
			return false;
		if (coverage == null) {
			if (other.coverage != null)
				return false;
		} else if (!coverage.equals(other.coverage))
			return false;
		if (gapLength == null) {
			if (other.gapLength != null)
				return false;
		} else if (!gapLength.equals(other.gapLength))
			return false;
		if (identity == null) {
			if (other.identity != null)
				return false;
		} else if (!identity.equals(other.identity))
			return false;
		if (initialShift == null) {
			if (other.initialShift != null)
				return false;
		} else if (!initialShift.equals(other.initialShift))
			return false;
		if (nNonSelfAligned == null) {
			if (other.nNonSelfAligned != null)
				return false;
		} else if (!nNonSelfAligned.equals(other.nNonSelfAligned))
			return false;
		if (rmsd == null) {
			if (other.rmsd != null)
				return false;
		} else if (!rmsd.equals(other.rmsd))
			return false;
		if (similarity == null) {
			if (other.similarity != null)
				return false;
		} else if (!similarity.equals(other.similarity))
			return false;
		if (tmScore == null) {
			if (other.tmScore != null)
				return false;
		} else if (!tmScore.equals(other.tmScore))
			return false;
		if (tmpr == null) {
			if (other.tmpr != null)
				return false;
		} else if (!tmpr.equals(other.tmpr))
			return false;
		if (zScore == null) {
			if (other.zScore != null)
				return false;
		} else if (!zScore.equals(other.zScore))
			return false;
		return true;
	}

	public Integer getAlignLength() {
		return alignLength;
	}

	public Float getAlignScore() {
		return alignScore;
	}

	public Float getAlternateTm() {
		return alternateTm;
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

	public Float getTmpr() {
		return tmpr;
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
		result = prime * result + (alternateTm == null ? 0 : alternateTm.hashCode());
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
		result = prime * result + (tmpr == null ? 0 : tmpr.hashCode());
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
		this.alternateTm = alternateTm;
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

	public void setTmpr(Float tmpr) {
		this.tmpr = tmpr;
	}

	public void setTmScore(Float tmScore) {
		this.tmScore = tmScore;
	}

	public void setzScore(Float zScore) {
		this.zScore = zScore;
	}

	@Override
	public String toString() {
		return "Alignment [alignLength=" + alignLength + ", alignScore=" + alignScore + ", block1Length="
				+ block1Length + ", block2Length=" + block2Length + ", coverage=" + coverage + ", gapLength="
				+ gapLength + ", identity=" + identity + ", rmsd=" + rmsd + ", similarity=" + similarity + ", tmScore="
				+ tmScore + ", alternateTm=" + alternateTm + ", zScore=" + zScore + ", initialShift=" + initialShift
				+ ", nNonSelfAligned=" + nNonSelfAligned + ", tmpr=" + tmpr + "]";
	}

}
