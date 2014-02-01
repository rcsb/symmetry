package org.biojava3.structure.align.symm.census2.analysis.ligands;

import java.io.Serializable;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import org.biojava.bio.structure.Atom;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

/**
 * A ligand that somehow relates to a symmetric structure.
 * 
 * @author dmyerstu
 */
public class Ligand implements Serializable {

	private static final long serialVersionUID = -3803769618055243227L;

	private double distanceToCentroid;
	private double distanceToAxis;
	private String formula;
	private boolean isMetallic;

	public Ligand(List<Atom> atoms, double distance, double distanceToAxis) {
		SortedMap<String, Integer> letters = new TreeMap<String, Integer>();
		for (Atom atom : atoms) {
			String name = atom.getElement().name();
			if (atom.getElement().getElementType().isMetal())
				isMetallic = true;
			StatUtils.plus(letters, name);
		}
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String, Integer> entry : letters.entrySet()) {
			sb.append(entry.getKey());
			if (entry.getValue() > 1)
				sb.append(entry.getValue());
		}
		this.formula = sb.toString();
		this.distanceToCentroid = distance;
		this.distanceToAxis = distanceToAxis;
	}
	
	public Ligand() {
		super();
	}

	public Ligand(String formula, boolean isMetal, double distance) {
		super();
		this.formula = formula;
		this.isMetallic = isMetal;
		this.distanceToCentroid = distance;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(distanceToAxis);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(distanceToCentroid);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((formula == null) ? 0 : formula.hashCode());
		result = prime * result + (isMetallic ? 1231 : 1237);
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
		Ligand other = (Ligand) obj;
		if (Double.doubleToLongBits(distanceToAxis) != Double
				.doubleToLongBits(other.distanceToAxis))
			return false;
		if (Double.doubleToLongBits(distanceToCentroid) != Double
				.doubleToLongBits(other.distanceToCentroid))
			return false;
		if (formula == null) {
			if (other.formula != null)
				return false;
		} else if (!formula.equals(other.formula))
			return false;
		if (isMetallic != other.isMetallic)
			return false;
		return true;
	}

	public String getFormula() {
		return formula;
	}

	public boolean isMetallic() {
		return isMetallic;
	}

	public void setFormula(String formula) {
		this.formula = formula;
	}

	public void setMetallic(boolean isMetallic) {
		this.isMetallic = isMetallic;
	}

	public double getDistanceToCentroid() {
		return distanceToCentroid;
	}

	public void setDistanceToCentroid(double distanceToCentroid) {
		this.distanceToCentroid = distanceToCentroid;
	}

	public double getDistanceToAxis() {
		return distanceToAxis;
	}

	public void setDistanceToAxis(double distanceToAxis) {
		this.distanceToAxis = distanceToAxis;
	}

	@Override
	public String toString() {
		return (isMetallic() ? "(metallic)" : "") + formula + " (" + StatUtils.formatD(distanceToCentroid) + "Å" + "to centroid;" + StatUtils.formatD(distanceToAxis) + "Å" + "to axis)";
	}

}
