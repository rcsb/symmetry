package org.biojava3.structure.quaternary.misc;

import java.util.Collections;
import java.util.List;

public class PdbBlastHit {
	private String pdbId = "";
    private List<String> chainIds = Collections.emptyList();
    private double eValue = 0;
    private double identity = 0;
    private int alignmentLength = 0;
    private int hitLength = 0;
    
    /**
	 * @return the pdbId
	 */
	public String getPdbId() {
		return pdbId;
	}
	/**
	 * @param pdbId the pdbId to set
	 */
	public void setPdbId(String pdbId) {
		this.pdbId = pdbId;
	}
	/**
	 * @return the chainIds
	 */
	public List<String> getChainIds() {
		return chainIds;
	}
	/**
	 * @param chainIds the chainIds to set
	 */
	public void setChainIds(List<String> chainIds) {
		this.chainIds = chainIds;
	}
	/**
	 * @return the eValue
	 */
	public double getEvalue() {
		return eValue;
	}
	/**
	 * @param eValue the eValue to set
	 */
	public void setEvalue(double eValue) {
		this.eValue = eValue;
	}
	/**
	 * @return the identity
	 */
	public double getIdentity() {
		return identity;
	}
	/**
	 * @param identity the identity to set
	 */
	public void setIdentity(double identity) {
		this.identity = identity;
	}  
	
	/**
	 * @return the alignmentLength
	 */
	public int getAlignmentLength() {
		return alignmentLength;
	}
	/**
	 * @param alignmentLength the alignmentLength to set
	 */
	public void setAlignmentLength(int alignmentLength) {
		this.alignmentLength = alignmentLength;
	}
	/**
	 * @return the hitLength
	 */
	public int getHitLength() {
		return hitLength;
	}
	/**
	 * @param hitLength the hitLength to set
	 */
	public void setHitLength(int hitLength) {
		this.hitLength = hitLength;
	}
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("PdbId: ");
		sb.append(pdbId);
		sb.append(" chain ids: ");
		sb.append(chainIds);
		sb.append(" eValue: ");
		sb.append(eValue);
		sb.append(" identity: ");
		sb.append(identity);
		sb.append(" alignment length: ");
		sb.append(alignmentLength);
		sb.append(" hit length: ");
		sb.append(hitLength);
		return sb.toString();
	}
}
