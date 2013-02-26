package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.Serializable;

import org.biojava3.structure.align.symm.census2.Alignment;
import org.biojava3.structure.align.symm.census2.Axis;
import org.biojava3.structure.align.symm.census2.Result;

public class Case implements Serializable {

	private static final long serialVersionUID = 1604599374800984540L;
	
	private int knownOrder;
	private String knownGroup;
	private Result result;
	public int getKnownOrder() {
		return knownOrder;
	}
	public void setKnownOrder(int knownOrder) {
		this.knownOrder = knownOrder;
	}
	public String getKnownGroup() {
		return knownGroup;
	}
	public void setKnownGroup(String knownGroup) {
		this.knownGroup = knownGroup;
	}
	public Result getResult() {
		return result;
	}
	public void setResult(Result result) {
		this.result = result;
	}
	@Override
	public String toString() {
		return "Case [knownOrder=" + knownOrder + ", knownGroup=" + knownGroup + ", result=" + result + "]";
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((knownGroup == null) ? 0 : knownGroup.hashCode());
		result = prime * result + knownOrder;
		result = prime * result + ((this.result == null) ? 0 : this.result.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		Case other = (Case) obj;
		if (knownGroup == null) {
			if (other.knownGroup != null) return false;
		} else if (!knownGroup.equals(other.knownGroup)) return false;
		if (knownOrder != other.knownOrder) return false;
		if (result == null) {
			if (other.result != null) return false;
		} else if (!result.equals(other.result)) return false;
		return true;
	}
	public Case(int knownOrder, String knownGroup, Result result) {
		this.knownOrder = knownOrder;
		this.knownGroup = knownGroup;
		this.result = result;
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
		return knownOrder > 1;
	}
	
}
