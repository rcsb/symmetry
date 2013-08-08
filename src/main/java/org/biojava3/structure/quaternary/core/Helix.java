
package org.biojava3.structure.quaternary.core;

import java.util.List;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;

/**
 *
 * @author Peter
 */
public class Helix {
    private double subunitRmsd = Double.MAX_VALUE;
    private double traceRmsd = Double.MAX_VALUE;
    private List<Integer> permutation;
    private List<List<Integer>> repeatUnits;
    private Matrix4d transformation;
    private double rise;
    private int nStart;
    private int fold;
    private int contacts;

    /**
     * @return the subunitRmsd
     */
    public double getSubunitRmsd() {
        return subunitRmsd;
    }

    /**
     * @param subunitRmsd the subunitRmsd to set
     */
    public void setSubunitRmsd(double subunitRmsd) {
        this.subunitRmsd = subunitRmsd;
    }

    /**
     * @return the traceRmsd
     */
    public double getTraceRmsd() {
        return traceRmsd;
    }

    /**
     * @param traceRmsd the traceRmsd to set
     */
    public void setTraceRmsd(double traceRmsd) {
        this.traceRmsd = traceRmsd;
    }

    /**
     * @return the permutation
     */
    public List<Integer> getPermutation() {
        return permutation;
    }

    /**
     * @param permutation the permutation to set
     */
    public void setPermutation(List<Integer> permutation) {
        this.permutation = permutation;
    }

    public List<List<Integer>> getRepeatUnits() {
		return repeatUnits;
	}

	public void setRepeatUnits(List<List<Integer>> repeatUnits) {
		this.repeatUnits = repeatUnits;
	}

	/**
     * @return the transformation
     */
    public Matrix4d getTransformation() {
        return transformation;
    }

    /**
     * @param transformation the transformation to set
     */
    public void setTransformation(Matrix4d transformation) {
        this.transformation = transformation;
    }

    public double getRise() {
		return rise;
	}

	public void setRise(double rise) {
		this.rise = rise;
	}
	
	/**
	 * Returns the pitch angle of the helix
	 * @param transformation helix transformation
	 * @return
	 */
	public double getAngle() {
		return getAxisAngle().angle;
	}
	
	/**
	 * Returns the AxisAngle of the helix transformation
	 * @param transformation helix transformation
	 * @return
	 */
	public AxisAngle4d getAxisAngle() {
		AxisAngle4d axis = new AxisAngle4d();
		axis.set(this.transformation);
		return axis;
	}

	public int getnStart() {
		return nStart;
	}

	public void setnStart(int nStart) {
		this.nStart = nStart;
	}

	/**
     * @return the fold
     */
    public int getFold() {
        return fold;
    }

    /**
     * @param fold the fold to set
     */
    public void setFold(int fold) {
        this.fold = fold;
    }
    
    public int getContacts() {
		return contacts;
	}

	public void setContacts(int contacts) {
		this.contacts = contacts;
	}

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Subunit RMSD  : " + getSubunitRmsd() + "\n");
        sb.append("CA RMSD       : " + getTraceRmsd() + "\n");
        sb.append("Permutation   : " + getPermutation() + "\n");
        sb.append("Repeat units  : " + getRepeatUnits() + "\n");
        sb.append("Rise          : " + getRise() + "\n");
        sb.append("Angle         : " + Math.toDegrees(getAngle()) +"\n");
        sb.append("Fold          : " + getFold() + "\n");
        return sb.toString();
    }
}
