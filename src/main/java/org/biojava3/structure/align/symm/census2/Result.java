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
package org.biojava3.structure.align.symm.census2;

import java.io.Serializable;

/**
 * a bean storing the results of the latest SCOP census.
 * 
 * @author Andreas Prlic
 * 
 */
public class Result implements Serializable {

	private static final long serialVersionUID = 2282745910618366982L;

	private Axis axis;
	private Alignment alignment;
	private String classification;
	private String description;
	private Boolean isSignificant;
	private Integer order;
	private String protodomain;
	private Integer rank;
	private String scopId;
	private Integer sunId;
	private Float fractionHelical;
	public Axis getAxis() {
		return axis;
	}
//	public Float getFractionHelical() {
//		return fractionHelical;
//	}
//	public void setFractionHelical(Float fractionHelical) {
//		this.fractionHelical = fractionHelical;
//	}
	public void setAxis(Axis axis) {
		this.axis = axis;
	}
	public Alignment getAlignment() {
		return alignment;
	}
	public void setAlignment(Alignment alignment) {
		this.alignment = alignment;
	}
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
	@Override
	public String toString() {
		return "Result [axis=" + axis + ", alignment=" + alignment + ", classification=" + classification
				+ ", description=" + description + ", isSignificant=" + isSignificant + ", order=" + order
				+ ", protodomain=" + protodomain + ", rank=" + rank + ", scopId=" + scopId + ", sunId=" + sunId
				+ ", fractionHelical=" + fractionHelical + "]";
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((alignment == null) ? 0 : alignment.hashCode());
		result = prime * result + ((axis == null) ? 0 : axis.hashCode());
		result = prime * result + ((classification == null) ? 0 : classification.hashCode());
		result = prime * result + ((description == null) ? 0 : description.hashCode());
		result = prime * result + ((fractionHelical == null) ? 0 : fractionHelical.hashCode());
		result = prime * result + ((isSignificant == null) ? 0 : isSignificant.hashCode());
		result = prime * result + ((order == null) ? 0 : order.hashCode());
		result = prime * result + ((protodomain == null) ? 0 : protodomain.hashCode());
		result = prime * result + ((rank == null) ? 0 : rank.hashCode());
		result = prime * result + ((scopId == null) ? 0 : scopId.hashCode());
		result = prime * result + ((sunId == null) ? 0 : sunId.hashCode());
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
		Result other = (Result) obj;
		if (alignment == null) {
			if (other.alignment != null)
				return false;
		} else if (!alignment.equals(other.alignment))
			return false;
		if (axis == null) {
			if (other.axis != null)
				return false;
		} else if (!axis.equals(other.axis))
			return false;
		if (classification == null) {
			if (other.classification != null)
				return false;
		} else if (!classification.equals(other.classification))
			return false;
		if (description == null) {
			if (other.description != null)
				return false;
		} else if (!description.equals(other.description))
			return false;
		if (fractionHelical == null) {
			if (other.fractionHelical != null)
				return false;
		} else if (!fractionHelical.equals(other.fractionHelical))
			return false;
		if (isSignificant == null) {
			if (other.isSignificant != null)
				return false;
		} else if (!isSignificant.equals(other.isSignificant))
			return false;
		if (order == null) {
			if (other.order != null)
				return false;
		} else if (!order.equals(other.order))
			return false;
		if (protodomain == null) {
			if (other.protodomain != null)
				return false;
		} else if (!protodomain.equals(other.protodomain))
			return false;
		if (rank == null) {
			if (other.rank != null)
				return false;
		} else if (!rank.equals(other.rank))
			return false;
		if (scopId == null) {
			if (other.scopId != null)
				return false;
		} else if (!scopId.equals(other.scopId))
			return false;
		if (sunId == null) {
			if (other.sunId != null)
				return false;
		} else if (!sunId.equals(other.sunId))
			return false;
		return true;
	}
	
}
