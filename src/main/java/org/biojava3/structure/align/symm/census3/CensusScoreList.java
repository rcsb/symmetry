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
package org.biojava3.structure.align.symm.census3;

import java.io.Serializable;

import javax.xml.bind.annotation.XmlTransient;

import org.biojava.bio.structure.align.model.AFPChain;

/**
 * The results of a symmetry alignment. Includes information for CE-Symm or
 * SymD.
 * 
 * @author dmyerstu
 */
public class CensusScoreList implements Serializable {

	private static final long serialVersionUID = 5560951230435020821L;

	private AdditionalScoreList additionalScoreList;

	private Integer alignLength;

	private Integer gapLength;

	private Float identity;

	private Float rmsd;

	private Float similarity;

	private Float tmScore;

	private Float zScore;

	public CensusScoreList() {

	}

	public CensusScoreList(AFPChain afpChain) {
		if (afpChain != null) {
			identity = (float) afpChain.getIdentity();
			similarity = (float) afpChain.getSimilarity();
			zScore = (float) afpChain.getProbability();
			rmsd = (float) afpChain.getTotalRmsdOpt();
			tmScore = (float) afpChain.getTMScore();
			alignLength = afpChain.getAlnLength();
			gapLength = afpChain.getGapLen();
		}
	}

	public CensusScoreList(Integer alignLength, Integer coverage, Float alignScore, Integer gapLength, Float identity,
			Float rmsd, Float similarity, Float tmScore, Float zScore, Integer block1Length, Integer block2Length) {
		this.alignLength = alignLength;
		this.gapLength = gapLength;
		this.identity = identity;
		this.rmsd = rmsd;
		this.similarity = similarity;
		this.tmScore = tmScore;
		this.zScore = zScore;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		CensusScoreList other = (CensusScoreList) obj;
		if (alignLength == null) {
			if (other.alignLength != null) return false;
		} else if (!alignLength.equals(other.alignLength)) return false;
		if (gapLength == null) {
			if (other.gapLength != null) return false;
		} else if (!gapLength.equals(other.gapLength)) return false;
		if (identity == null) {
			if (other.identity != null) return false;
		} else if (!identity.equals(other.identity)) return false;
		if (rmsd == null) {
			if (other.rmsd != null) return false;
		} else if (!rmsd.equals(other.rmsd)) return false;
		if (similarity == null) {
			if (other.similarity != null) return false;
		} else if (!similarity.equals(other.similarity)) return false;
		if (tmScore == null) {
			if (other.tmScore != null) return false;
		} else if (!tmScore.equals(other.tmScore)) return false;
		if (zScore == null) {
			if (other.zScore != null) return false;
		} else if (!zScore.equals(other.zScore)) return false;
		return true;
	}

	// For now, ignore extra scores in XML since it's abstract.
	@XmlTransient
	// Could also list all possible concrete implementations here to include it.
	//@XmlElements({
	//	@XmlElement(type=SinglePropScoreList.class),
	//	@XmlElement(type=MapScoreList.class),
	//})
	public AdditionalScoreList getAdditionalScoreList() {
		return additionalScoreList;
	}

	public Integer getAlignLength() {
		return alignLength;
	}

	public Integer getGapLength() {
		return gapLength;
	}

	public Float getIdentity() {
		return identity;
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
		result = prime * result + (gapLength == null ? 0 : gapLength.hashCode());
		result = prime * result + (identity == null ? 0 : identity.hashCode());
		result = prime * result + (rmsd == null ? 0 : rmsd.hashCode());
		result = prime * result + (similarity == null ? 0 : similarity.hashCode());
		result = prime * result + (tmScore == null ? 0 : tmScore.hashCode());
		result = prime * result + (zScore == null ? 0 : zScore.hashCode());
		return result;
	}

	public void setAdditionalScoreList(AdditionalScoreList additionalScoreList) {
		this.additionalScoreList = additionalScoreList;
	}

	public void setAlignLength(Integer alignLength) {
		this.alignLength = alignLength;
	}

	public void setGapLength(Integer gapLength) {
		this.gapLength = gapLength;
	}

	public void setIdentity(Float identity) {
		this.identity = identity;
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

	@Override
	public String toString() {
		return "CensusScoreList [alignLength=" + alignLength + ", gapLength=" + gapLength + ", identity=" + identity
				+ ", rmsd=" + rmsd + ", similarity=" + similarity + ", tmScore=" + tmScore + ", zScore=" + zScore + "]";
	}

}
