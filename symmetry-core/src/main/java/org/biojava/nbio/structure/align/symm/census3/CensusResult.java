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
package org.biojava.nbio.structure.align.symm.census3;

import java.io.Serializable;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;

/**
 * An entry in a symmetry census.
 * @author dmyersturnbull
 */
public class CensusResult implements Serializable {

	private static final long serialVersionUID = 2282745910618366982L;

	private CensusAlignment alignment;
	private CensusAxis axis;
	private CensusSymmetryGroup group;
	private String id;
	
	private CensusScoreList scoreList;
	private String alignedUnit;

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		CensusResult other = (CensusResult) obj;
		if (scoreList == null) {
			if (other.scoreList != null) return false;
		} else if (!scoreList.equals(other.scoreList)) return false;
		if (axis == null) {
			if (other.axis != null) return false;
		} else if (!axis.equals(other.axis)) return false;
		if (alignedUnit == null) {
			if (other.alignedUnit != null) return false;
		} else if (!alignedUnit.equals(other.alignedUnit)) return false;
		if (id == null) {
			if (other.id != null) return false;
		} else if (!id.equals(other.id)) return false;
		return true;
	}

	@XmlElement(name = "scores")
	public CensusScoreList getScoreList() {
		return scoreList;
	}

	public CensusAlignment getAlignment() {
		return alignment;
	}

	public CensusAxis getAxis() {
		return axis;
	}

	public CensusSymmetryGroup getGroup() {
		return group;
	}

	public String getId() {
		return id;
	}

	public String getAlignedUnit() {
		return alignedUnit;
	}

	public int getOrder() {
		if (group == null) return 1;
		return group.getOrder();
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (scoreList == null ? 0 : scoreList.hashCode());
		result = prime * result + (axis == null ? 0 : axis.hashCode());
		result = prime * result + (group == null ? 0 : group.hashCode());
		result = prime * result + (alignedUnit == null ? 0 : alignedUnit.hashCode());
		result = prime * result + (id == null ? 0 : id.hashCode());
		return result;
	}

	public void setScoreList(CensusScoreList scoreList) {
		this.scoreList = scoreList;
	}

	public void setAlignment(CensusAlignment alignment) {
		this.alignment = alignment;
	}

	public void setAxis(CensusAxis axis) {
		this.axis = axis;
	}

	public void setGroup(CensusSymmetryGroup group) {
		this.group = group;
	}

	@XmlAttribute
	public void setId(String id) {
		this.id = id;
	}

	public void setAlignedUnit(String alignedUnit) {
		this.alignedUnit = alignedUnit;
	}

	@Override
	public String toString() {
		return "CensusResult [scoreList=" + scoreList + ", alignment=" + alignment + ", axis=" + axis + ", group="
				+ group + ", alignedUnit=" + alignedUnit + ", id=" + id + "]";
	}

}
