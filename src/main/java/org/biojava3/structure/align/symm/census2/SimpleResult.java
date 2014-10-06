/**
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
 * Created on 25 Mar 2013
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava3.structure.align.symm.census2;

/**
 * For now, just used to represent a {@link Result} as a simple bean with only primitives and strings as fields.
 * This was needed for the servlet.
 * @author dmyersturnbull
 * @deprecated
 */
@Deprecated
public class SimpleResult {

	public SimpleResult(Result r) {
		classification = r.getClassification();
		description = r.getDescription();
		isSignificant = r.getIsSignificant();
		order = r.getOrder();
		protodomain = r.getProtodomain();
		rank = r.getRank();
		scopId = r.getScopId();
		sunId = r.getSunId();
		if (r.getAxis() != null) {
			theta = r.getAxis().getTheta();
		}
		if (r.getAlignment() != null) {
			tmScore = r.getAlignment().getTmScore();
			zScore = r.getAlignment().getzScore();
			alignScore = r.getAlignment().getAlignScore();
			gapLength = r.getAlignment().getGapLength();
			identity = r.getAlignment().getIdentity();
			similarity = r.getAlignment().getSimilarity();
			alignLength = r.getAlignment().getAlignLength();
			rmsd = r.getAlignment().getRmsd();
			length1 = gapLength + alignLength;
		}
	}

	public SimpleResult() {
		super();
	}

	private String classification;
	private String description;
	private Boolean isSignificant;
	private Integer order;
	private String protodomain;
	private Integer rank;
	private String scopId;
	private Integer sunId;
	private Float theta;
	private Float tmScore;
	private Float zScore;
	private Float alignScore;
	private Integer gapLength;
	private Float identity;
	private Float similarity;
	private Integer alignLength;
	private Float rmsd;
	private Integer length1;

	public String getClassification() {
		return classification;
	}
	public void setClassification(String classification) {
		this.classification = classification;
	}
	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
	}
	public Boolean getIsSignificant() {
		return isSignificant;
	}
	public void setIsSignificant(Boolean isSignificant) {
		this.isSignificant = isSignificant;
	}
	public Integer getOrder() {
		return order;
	}
	public void setOrder(Integer order) {
		this.order = order;
	}
	public String getProtodomain() {
		return protodomain;
	}
	public void setProtodomain(String protodomain) {
		this.protodomain = protodomain;
	}
	public Integer getRank() {
		return rank;
	}
	public void setRank(Integer rank) {
		this.rank = rank;
	}
	public String getScopId() {
		return scopId;
	}
	public void setScopId(String scopId) {
		this.scopId = scopId;
	}
	public Integer getSunId() {
		return sunId;
	}
	public void setSunId(Integer sunId) {
		this.sunId = sunId;
	}
	public Float getTheta() {
		return theta;
	}
	public void setTheta(Float theta) {
		this.theta = theta;
	}
	public Float getTmScore() {
		return tmScore;
	}
	public void setTmScore(Float tmScore) {
		this.tmScore = tmScore;
	}
	public Float getzScore() {
		return zScore;
	}
	public void setzScore(Float zScore) {
		this.zScore = zScore;
	}
	public Float getAlignScore() {
		return alignScore;
	}
	public void setAlignScore(Float alignScore) {
		this.alignScore = alignScore;
	}
	public Integer getGapLength() {
		return gapLength;
	}
	public void setGapLength(Integer gapLength) {
		this.gapLength = gapLength;
	}
	public Float getIdentity() {
		return identity;
	}
	public void setIdentity(Float identity) {
		this.identity = identity;
	}
	public Float getSimilarity() {
		return similarity;
	}
	public void setSimilarity(Float similarity) {
		this.similarity = similarity;
	}
	public Integer getAlignLength() {
		return alignLength;
	}
	public void setAlignLength(Integer alignLength) {
		this.alignLength = alignLength;
	}
	public Float getRmsd() {
		return rmsd;
	}
	public void setRmsd(Float rmsd) {
		this.rmsd = rmsd;
	}
	public Integer getLength1() {
		return length1;
	}
	public void setLength1(Integer length1) {
		this.length1 = length1;
	}

}
