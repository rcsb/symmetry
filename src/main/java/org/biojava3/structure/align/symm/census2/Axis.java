/*
 * BioJava development code
 * 
 * This code may be freely distributed and modified under the terms of the GNU Lesser General Public Licence. This
 * should be distributed with the code. If you do not have a copy, see:
 * 
 * http://www.gnu.org/copyleft/lesser.html
 * 
 * Copyright for this code is held jointly by the individual authors. These should be listed in @author doc comments.
 * 
 * For more information on the BioJava project and its aims, or to join the biojava-l mailing list, visit the home page
 * at:
 * 
 * http://www.biojava.org/
 * 
 * Created on 2013-02-22
 */
package org.biojava3.structure.align.symm.census2;

import java.io.Serializable;

import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.align.util.RotationAxis;

/**
 * An axis of rotation.
 * 
 * @author dmyerstu
 */
public class Axis implements Serializable {

	private static final long serialVersionUID = -1523140268972093071L;

	private Float orthogonal;
	private Float screw;
	private Float theta;
	private Float x;
	private Float y;
	private Float z;

	public Axis() {
	}

	public Axis(Float angle, Float x, Float y, Float z, Float screw, Float orthogonal) {
		theta = angle;
		this.x = x;
		this.y = y;
		this.z = z;
		this.screw = screw;
		this.orthogonal = orthogonal;
	}

	public Axis(RotationAxis rot) {
		theta = (float) rot.getAngle();
		x = (float) rot.getRotationAxis().getX();
		y = (float) rot.getRotationAxis().getY();
		z = (float) rot.getRotationAxis().getZ();
		screw = (float) (Calc.amount(rot.getScrewTranslation()) / Calc.amount(rot.getRotationAxis()));
		orthogonal = (float) (Calc.amount(rot.getOtherTranslation()) / Calc.amount(rot.getRotationAxis()));
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		Axis other = (Axis) obj;
		if (Double.doubleToLongBits(theta) != Double.doubleToLongBits(other.theta)) return false;
		if (Double.doubleToLongBits(y) != Double.doubleToLongBits(other.y)) return false;
		if (Double.doubleToLongBits(z) != Double.doubleToLongBits(other.z)) return false;
		if (Double.doubleToLongBits(screw) != Double.doubleToLongBits(other.screw)) return false;
		if (Double.doubleToLongBits(orthogonal) != Double.doubleToLongBits(other.orthogonal)) return false;
		if (Double.doubleToLongBits(x) != Double.doubleToLongBits(other.x)) return false;
		return true;
	}

	public Double evaluateEpsilon(int order) {
		if (order < 2) return null;
		final double gamma = 2 * Math.PI / order;
		double theta = Math.abs(getTheta());
		double minEpsilon = Math.abs(theta - gamma);
		double delta = gamma;
		for (int k = 1; delta <= theta; k++, delta += gamma) {
			if (order % k != 0) continue;
			double epsilon = Math.abs(theta - delta);
			if (epsilon < minEpsilon) minEpsilon = epsilon;
		}
		return minEpsilon;
	}

	public Float getOrthogonal() {
		return orthogonal;
	}

	public Float getScrew() {
		return screw;
	}

	public Float getTheta() {
		return theta;
	}

	public Float getX() {
		return x;
	}

	public Float getY() {
		return y;
	}

	public Float getZ() {
		return z;
	}

	public int guessOrder() {
		return guessOrder(1.0 * Math.PI / 180, 8);
	}

	/**
	 * Guesses an order of rotational symmetry from {@link #getTheta() theta}.
	 * 
	 * @param threshold
	 *            The maximal deviation of theta from the expected angle (2π/order) required
	 * @param maxOrder
	 *            The maximum order of symmetry to consider
	 * @return The order between 2 and {@code maxOrder} that results in the lowest deviation |θ − 2π/order|, or
	 *         {@code 1} if no deviation less than {@code threshold} is found
	 */
	public int guessOrder(double threshold, int maxOrder) {
		double bestDelta = threshold;
		int bestOrder = 1;
		for (int order = 2; order < maxOrder; order++) {
			double delta = Math.abs(2 * Math.PI / order - theta);
			if (delta < bestDelta) {
				bestOrder = order;
				bestDelta = delta;
			}
		}
		return bestOrder;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(theta);
		result = prime * result + (int) (temp ^ temp >>> 32);
		temp = Double.doubleToLongBits(y);
		result = prime * result + (int) (temp ^ temp >>> 32);
		temp = Double.doubleToLongBits(z);
		result = prime * result + (int) (temp ^ temp >>> 32);
		temp = Double.doubleToLongBits(screw);
		result = prime * result + (int) (temp ^ temp >>> 32);
		temp = Double.doubleToLongBits(orthogonal);
		result = prime * result + (int) (temp ^ temp >>> 32);
		temp = Double.doubleToLongBits(x);
		result = prime * result + (int) (temp ^ temp >>> 32);
		return result;
	}

	public void setOrthogonal(Float orthogonal) {
		this.orthogonal = orthogonal;
	}

	public void setScrew(Float screw) {
		this.screw = screw;
	}

	public void setTheta(Float theta) {
		this.theta = theta;
	}

	public void setX(Float x) {
		this.x = x;
	}

	public void setY(Float y) {
		this.y = y;
	}

	public void setZ(Float z) {
		this.z = z;
	}

	@Override
	public String toString() {
		return "Axis [angle=" + theta + ", x=" + x + ", y=" + y + ", z=" + z + ", screw=" + screw + ", orthogonal="
				+ orthogonal + "]";
	}

}
