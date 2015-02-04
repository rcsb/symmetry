package org.biojava.nbio.ramachandran;

public class ScoreMap implements Comparable<ScoreMap>{
	Integer phi;
	Integer psi;
	int count;
	public ScoreMap() {
	
		count = 0;
	}
	public Integer getPhi() {
		return phi;
	}
	public void setPhi(Integer phi) {
		this.phi = phi;
	}
	public Integer getPsi() {
		return psi;
	}
	public void setPsi(Integer psi) {
		this.psi = psi;
	}
	public int getCount() {
		return count;
	}
	public void setCount(int count) {
		this.count = count;
	}
	@Override
	public int compareTo(ScoreMap o) {
		if ( this.equals(o))
			return 0;
		
		if ( o.getPhi().equals(phi)){		
			return o.getPsi().compareTo(psi);
		} else {
			return o.getPhi().compareTo(phi);
		}
	}
	
	
	public void incrementCount() {
		count++;
		
	}
	@Override
	public String toString() {
		return "ScoreMap [phi=" + phi + ", psi=" + psi + ", count=" + count
				+ "]";
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + count;
		result = prime * result + ((phi == null) ? 0 : phi.hashCode());
		result = prime * result + ((psi == null) ? 0 : psi.hashCode());
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
		ScoreMap other = (ScoreMap) obj;
		
		if (phi == null) {
			if (other.phi != null)
				return false;
		} else if (!phi.equals(other.phi)) {
			if (  ! overlap(phi,other.phi) )
				return false;
		}
		if (psi == null) {
			if (other.psi != null)
				return false;
		} else if (!psi.equals(other.psi)) {
			if (  ! overlap(psi,other.psi) )
			return false;
		}
		return true;
	}
	static final int w = 2;
	private boolean overlap(Integer a, Integer b) {
		
		if (a >= b- w && a <= b + w )
			return true;
		
		return false;
	}
	
	
	
}
