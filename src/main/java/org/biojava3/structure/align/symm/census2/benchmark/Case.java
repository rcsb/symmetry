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
package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.Serializable;

import org.biojava3.structure.align.symm.census2.Alignment;
import org.biojava3.structure.align.symm.census2.Axis;
import org.biojava3.structure.align.symm.census2.Result;

/**
 * A single case (result) to benchmark. Contains a {@link Result} and a known {@link #getKnownOrder() order} and {@link #getKnownGroup() group}.
 * @author dmyerstu
 *
 */
public class Case implements Serializable {

	private static final long serialVersionUID = 1604599374800984540L;
	
	private KnownInfo knownInfo;
	private Result result;
	public int getKnownOrder() {
		return knownInfo.getOrder();
	}
	public String getKnownGroup() {
		return knownInfo.getGroup();
	}
	
	public KnownInfo getKnownInfo() {
		return knownInfo;
	}
	public void setKnownInfo(KnownInfo knownInfo) {
		this.knownInfo = knownInfo;
	}
	public Result getResult() {
		return result;
	}
	public void setResult(Result result) {
		this.result = result;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((knownInfo == null) ? 0 : knownInfo.hashCode());
		result = prime * result + ((this.result == null) ? 0 : this.result.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		Case other = (Case) obj;
		if (knownInfo == null) {
			if (other.knownInfo != null) return false;
		} else if (!knownInfo.equals(other.knownInfo)) return false;
		if (result == null) {
			if (other.result != null) return false;
		} else if (!result.equals(other.result)) return false;
		return true;
	}
	
	@Override
	public String toString() {
		return "Case [knownInfo=" + knownInfo + ", result=" + result + "]";
	}
	public Case() {
	}
	public Axis getAxis() {
		return result.getAxis();
	}
	public Alignment getAlignment() {
		return result.getAlignment();
	}
	public String getClassification() {
		return result.getClassification();
	}
	public String getDescription() {
		return result.getDescription();
	}
	public Integer getOrder() {
		return result.getOrder();
	}
	public String getProtodomain() {
		return result.getProtodomain();
	}
	public Integer getRank() {
		return result.getRank();
	}
	public String getScopId() {
		return result.getScopId();
	}
	public Integer getSunId() {
		return result.getSunId();
	}
	public boolean hasKnownSymmetry() {
		return getKnownOrder() > 1;
	}
	public boolean hasDecidedSymmetry() {
		return getOrder() > 1;
	}
}
