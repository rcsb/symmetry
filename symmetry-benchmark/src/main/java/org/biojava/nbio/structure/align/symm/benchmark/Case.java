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
package org.biojava.nbio.structure.align.symm.benchmark;

import java.io.Serializable;

import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.symm.census3.CensusAxis;
import org.biojava.nbio.structure.align.symm.census3.CensusResult;
import org.biojava.nbio.structure.align.symm.census3.CensusScoreList;


/**
 * A single case (result) to benchmark. Contains a {@link Result} (see symmetry project) and a known space group as {@link KnownInfo}.
 * @author dmyerstu
 * @see Sample, which contains many Cases
 */
public class Case implements Serializable {

	private static final long serialVersionUID = 1604599374800984540L;
	
	private KnownInfo knownInfo;
	private CensusResult result;
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
	public CensusResult getResult() {
		return result;
	}
	public void setResult(CensusResult result) {
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
	public CensusAxis getAxis() {
		return result.getAxis();
	}
	public CensusScoreList getScoreList() {
		return result.getScoreList();
	}
	public Integer getOrder() {
		return result.getOrder();
	}
	public String getAlignedUnit() {
		return result.getAlignedUnit();
	}
	public String getScopId() {
		return result.getId();
	}
}
