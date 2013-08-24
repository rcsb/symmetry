/**
 * 
 */
package org.biojava3.structure.quaternary.core;

/**
 * @author Peter
 *
 */
public class QuatSymmetryScores {
	private double minRmsd = 0;
	private double maxRmsd = 0;
	private double rmsd = 0;
	private double minTm = 1;
	private double maxTm = 1;
	private double tm = 1;
	private double rmsdCenters = 0;
	
	/**
	 * @return the minRmsd
	 */
	public double getMinRmsd() {
		return minRmsd;
	}
	/**
	 * @param minRmsd the minRmsd to set
	 */
	public void setMinRmsd(double minRmsd) {
		this.minRmsd = minRmsd;
	}
	/**
	 * @return the maxRmsd
	 */
	public double getMaxRmsd() {
		return maxRmsd;
	}
	/**
	 * @param maxRmsd the maxRmsd to set
	 */
	public void setMaxRmsd(double maxRmsd) {
		this.maxRmsd = maxRmsd;
	}
	/**
	 * @return the rmsd
	 */
	public double getRmsd() {
		return rmsd;
	}
	/**
	 * @param rmsd the rmsd to set
	 */
	public void setRmsd(double rmsd) {
		this.rmsd = rmsd;
	}
	/**
	 * @return the minTm
	 */
	public double getMinTm() {
		return minTm;
	}
	/**
	 * @param minTm the minTm to set
	 */
	public void setMinTm(double minTm) {
		this.minTm = minTm;
	}
	/**
	 * @return the maxTm
	 */
	public double getMaxTm() {
		return maxTm;
	}
	/**
	 * @param maxTm the maxTm to set
	 */
	public void setMaxTm(double maxTm) {
		this.maxTm = maxTm;
	}
	/**
	 * @return the tm
	 */
	public double getTm() {
		return tm;
	}
	/**
	 * @param tm the tm to set
	 */
	public void setTm(double tm) {
		this.tm = tm;
	}
	/**
	 * @return the rmsdCenters
	 */
	public double getRmsdCenters() {
		return rmsdCenters;
	}
	/**
	 * @param rmsdCenters the rmsdCenters to set
	 */
	public void setRmsdCenters(double rmsdCenters) {
		this.rmsdCenters = rmsdCenters;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("minimum RMSD: ");
		sb.append(getMinRmsd());
		sb.append("\n");
		sb.append("maximum RMSD: ");
		sb.append(getMaxRmsd());
		sb.append("\n");
		sb.append("RMSD        : ");
		sb.append(getRmsd());
		sb.append("\n");
		sb.append("minimum TM  : ");
		sb.append(getMinTm());
		sb.append("\n");
		sb.append("maximum TM  : ");
		sb.append(getMaxTm());
		sb.append("\n");
		sb.append("TM          : ");
		sb.append(getTm());
		sb.append("\n");
		sb.append("center RMSD: ");
		sb.append(getRmsdCenters());
		sb.append("\n");

		return sb.toString();
	}
}
