package org.biojava3.structure.align.symm.census2.analysis;

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

	private double distance;
	private String formula;
	private boolean isMetallic;

	public Ligand(List<Atom> atoms, double distance) {
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
		this.distance = distance;
	}
	
	public Ligand() {
		super();
	}

	public Ligand(String formula, boolean isMetal, double distance) {
		super();
		this.formula = formula;
		this.isMetallic = isMetal;
		this.distance = distance;
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
		if (Double.doubleToLongBits(distance) != Double.doubleToLongBits(other.distance))
			return false;
		if (formula == null) {
			if (other.formula != null)
				return false;
		} else if (!formula.equals(other.formula))
			return false;
		return true;
	}

	public double getDistance() {
		return distance;
	}

	public String getFormula() {
		return formula;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(distance);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((formula == null) ? 0 : formula.hashCode());
		return result;
	}

	public boolean isMetallic() {
		return isMetallic;
	}

	public void setDistance(double distance) {
		this.distance = distance;
	}

	public void setFormula(String formula) {
		this.formula = formula;
	}

	public void setMetallic(boolean isMetallic) {
		this.isMetallic = isMetallic;
	}

	@Override
	public String toString() {
		return (isMetallic() ? "*" : "") + formula + " (" + StatUtils.formatD(distance) + "Å" + ")";
	}

}
