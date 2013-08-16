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
 * Created on Sep 30, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava3.structure.align.symm.census;

import java.io.File;
import java.io.Serializable;

import org.biojava.bio.structure.align.util.AtomCache;

/** a bean storing the results of the latest SCOP census.
 * 
 * @author Andreas Prlic
 * 
 * @deprecated Replaced by {@link Results}
 *
 */
@Deprecated
public class CensusResult implements Serializable  {

	/**
	 * 
	 */
	private static final long serialVersionUID = 2282745910618366982L;


	Integer rank;
	Boolean isSignificant;
	String classificationId;
	String name;
	Double zScore;
	Double rmsd;
	Double tmScore;
	Double alignScore;
	Float identity;
	Float similarity;
	Integer length1;
	Integer aligLength;
	Float angle;
	Integer order;
	String description;
	Character scopClass;
	String protoDomain;
	Integer sunid;
	
	public Integer getRank() {
		return rank;
	}
	public void setRank(Integer rank) {
		this.rank = rank;
	}
	public Boolean getIsSignificant() {
		return isSignificant;
	}
	public void setIsSignificant(Boolean isSignificant) {
		this.isSignificant = isSignificant;
	}
	public String getClassificationId() {
		return classificationId;
	}
	public void setClassificationId(String classificationId) {
		this.classificationId = classificationId;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public Double getzScore() {
		return zScore;
	}
	public void setzScore(Double zScore) {
		this.zScore = zScore;
	}
	public Double getRmsd() {
		return rmsd;
	}
	public void setRmsd(Double rmsd) {
		this.rmsd = rmsd;
	}
	public Double getTmScore() {
		return tmScore;
	}
	public void setTmScore(Double tmScore) {
		this.tmScore = tmScore;
	}
	public Double getAlignScore() {
		return alignScore;
	}
	public void setAlignScore(Double alignScore) {
		this.alignScore = alignScore;
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
	public Integer getLength1() {
		return length1;
	}
	public void setLength1(Integer length1) {
		this.length1 = length1;
	}
	public Integer getAligLength() {
		return aligLength;
	}
	public void setAligLength(Integer aligLength) {
		this.aligLength = aligLength;
	}
	public Float getAngle() {
		return angle;
	}
	public void setAngle(Float angle) {
		this.angle = angle;
	}
	public Integer getOrder() {
		return order;
	}
	public void setOrder(Integer order) {
		this.order = order;
	}
	
	
	
	
	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
	}
	
	
	
	
	public Character getScopClass() {
		return scopClass;
	}
	public void setScopClass(Character scopClass) {
		this.scopClass = scopClass;
	}
	
	
	
	public String getProtoDomain() {
		return protoDomain;
	}
	public void setProtoDomain(String protoDomain) {
		this.protoDomain = protoDomain;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((aligLength == null) ? 0 : aligLength.hashCode());
		result = prime * result
				+ ((alignScore == null) ? 0 : alignScore.hashCode());
		result = prime * result + ((angle == null) ? 0 : angle.hashCode());
		result = prime
				* result
				+ ((classificationId == null) ? 0 : classificationId.hashCode());
		result = prime * result
				+ ((identity == null) ? 0 : identity.hashCode());
		result = prime * result
				+ ((isSignificant == null) ? 0 : isSignificant.hashCode());
		result = prime * result + ((length1 == null) ? 0 : length1.hashCode());
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result + ((order == null) ? 0 : order.hashCode());
		result = prime * result + ((rank == null) ? 0 : rank.hashCode());
		result = prime * result + ((rmsd == null) ? 0 : rmsd.hashCode());
		result = prime * result
				+ ((similarity == null) ? 0 : similarity.hashCode());
		result = prime * result + ((tmScore == null) ? 0 : tmScore.hashCode());
		result = prime * result + ((zScore == null) ? 0 : zScore.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		CensusResult other = (CensusResult) obj;
		if (aligLength == null) {
			if (other.aligLength != null)
				return false;
		} else if (!aligLength.equals(other.aligLength))
			return false;
		if (alignScore == null) {
			if (other.alignScore != null)
				return false;
		} else if (!alignScore.equals(other.alignScore))
			return false;
		if (angle == null) {
			if (other.angle != null)
				return false;
		} else if (!angle.equals(other.angle))
			return false;
		if (classificationId == null) {
			if (other.classificationId != null)
				return false;
		} else if (!classificationId.equals(other.classificationId))
			return false;
		if (identity == null) {
			if (other.identity != null)
				return false;
		} else if (!identity.equals(other.identity))
			return false;
		if (isSignificant == null) {
			if (other.isSignificant != null)
				return false;
		} else if (!isSignificant.equals(other.isSignificant))
			return false;
		if (length1 == null) {
			if (other.length1 != null)
				return false;
		} else if (!length1.equals(other.length1))
			return false;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (order == null) {
			if (other.order != null)
				return false;
		} else if (!order.equals(other.order))
			return false;
		if (rank == null) {
			if (other.rank != null)
				return false;
		} else if (!rank.equals(other.rank))
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
		if (zScore == null) {
			if (other.zScore != null)
				return false;
		} else if (!zScore.equals(other.zScore))
			return false;
		return true;
	}
	@Override
	public String toString() {
		return "CensusResult [rank=" + rank + ", isSignificant="
				+ isSignificant + ", classificationId=" + classificationId
				+ ", name=" + name + ", zScore=" + zScore + ", rmsd=" + rmsd
				+ ", tmScore=" + tmScore + ", alignScore=" + alignScore
				+ ", identity=" + identity + ", similarity=" + similarity
				+ ", length1=" + length1 + ", aligLength=" + aligLength
				+ ", angle=" + angle + ", order=" + order + ", description="
				+ description + ", scopClass=" + scopClass + ", getRank()="
				+ getRank() + ", getIsSignificant()=" + getIsSignificant()
				+ ", getClassificationId()=" + getClassificationId()
				+ ", getName()=" + getName() + ", getzScore()=" + getzScore()
				+ ", getRmsd()=" + getRmsd() + ", getTmScore()=" + getTmScore()
				+ ", getAlignScore()=" + getAlignScore() + ", getIdentity()="
				+ getIdentity() + ", getSimilarity()=" + getSimilarity()
				+ ", getLength1()=" + getLength1() + ", getAligLength()="
				+ getAligLength() + ", getAngle()=" + getAngle()
				+ ", getOrder()=" + getOrder() + ", getDescription()="
				+ getDescription() + ", getScopClass()=" + getScopClass() 
				+ getProtoDomain() 
				+ "]";
	}
	
	
	public boolean hasProtoDomainScanReady(AtomCache cache){
		
		if ( ! isSignificant )
			return false;
				
		String resultFile = CallableProtodomainComparison.getResultFilePath(cache, protoDomain);
		File f = new File(resultFile);
		return f.exists();
		
	}
	public Integer getSunid() {
		return sunid;
	}
	public void setSunid(Integer sunid) {
		this.sunid = sunid;
	}
	
	

}
