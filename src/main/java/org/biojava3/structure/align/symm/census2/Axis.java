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
package org.biojava3.structure.align.symm.census2;

import java.io.Serializable;

import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.gui.RotationAxis;

/**
 * An axis of rotation.
 * @author dmyerstu
 */
public class Axis implements Serializable {

	private static final long serialVersionUID = -1523140268972093071L;
	
	private Float theta;
	private Float x;
	private Float y;
	private Float z;
	private Float screw;
	private Float orthogonal;
	public Axis(Float angle, Float yaw, Float pitch, Float roll, Float screw, Float skew) {
		this.theta = angle;
		this.x = yaw;
		this.y = pitch;
		this.z = roll;
		this.screw = screw;
		this.orthogonal = skew;
	}
	
	public Float getTheta() {
		return theta;
	}

	public void setTheta(Float theta) {
		this.theta = theta;
	}

	public Float getX() {
		return x;
	}

	public void setX(Float x) {
		this.x = x;
	}

	public Float getY() {
		return y;
	}

	public void setY(Float y) {
		this.y = y;
	}

	public Float getZ() {
		return z;
	}

	public void setZ(Float z) {
		this.z = z;
	}

	public Float getScrew() {
		return screw;
	}

	public void setScrew(Float screw) {
		this.screw = screw;
	}

	public Float getOrthogonal() {
		return orthogonal;
	}

	public void setOrthogonal(Float orthogonal) {
		this.orthogonal = orthogonal;
	}

	public Axis() {
	}
	public Axis(RotationAxis rot) {
		theta = (float) rot.getTheta();
		x = (float) rot.getRotationAxis().getX();
		y = (float) rot.getRotationAxis().getY();
		z = (float) rot.getRotationAxis().getZ();
		screw = (float) (Calc.amount(rot.getScrewTranslation()) / Calc.amount(rot.getRotationAxis()));
		orthogonal = (float) (Calc.amount(rot.getOtherTranslation()) / Calc.amount(rot.getRotationAxis()));
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(theta);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(y);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(z);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(screw);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(orthogonal);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(x);
		result = prime * result + (int) (temp ^ (temp >>> 32));
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
		Axis other = (Axis) obj;
		if (Double.doubleToLongBits(theta) != Double.doubleToLongBits(other.theta))
			return false;
		if (Double.doubleToLongBits(y) != Double.doubleToLongBits(other.y))
			return false;
		if (Double.doubleToLongBits(z) != Double.doubleToLongBits(other.z))
			return false;
		if (Double.doubleToLongBits(screw) != Double.doubleToLongBits(other.screw))
			return false;
		if (Double.doubleToLongBits(orthogonal) != Double.doubleToLongBits(other.orthogonal))
			return false;
		if (Double.doubleToLongBits(x) != Double.doubleToLongBits(other.x))
			return false;
		return true;
	}
	@Override
	public String toString() {
		return "Axis [angle=" + theta + ", x=" + x + ", y=" + y + ", z=" + z + ", screw=" + screw
				+ ", orthogonal=" + orthogonal + "]";
	}
	
}
