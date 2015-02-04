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
package org.biojava.nbio.structure.align.symm.census3;

import java.io.Serializable;

import javax.xml.bind.annotation.XmlAttribute;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.align.util.RotationAxis;

/**
 * An axis of rotation.
 * 
 * @author dmyersturnbull
 */
public class CensusAxis implements Serializable {

	public static class SymmetryAxisVector {

		private Float x;

		private Float y;

		private Float z;

		public SymmetryAxisVector() {
			this(0f,0f,0f);
		}

		public SymmetryAxisVector(Float x, Float y, Float z) {
			super();
			this.x = x;
			this.y = y;
			this.z = z;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			SymmetryAxisVector other = (SymmetryAxisVector) obj;
			if (x == null) {
				if (other.x != null) return false;
			} else if (!x.equals(other.x)) return false;
			if (y == null) {
				if (other.y != null) return false;
			} else if (!y.equals(other.y)) return false;
			if (z == null) {
				if (other.z != null) return false;
			} else if (!z.equals(other.z)) return false;
			return true;
		}

		@XmlAttribute
		public Float getX() {
			return x;
		}

		@XmlAttribute
		public Float getY() {
			return y;
		}

		@XmlAttribute
		public Float getZ() {
			return z;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (x == null ? 0 : x.hashCode());
			result = prime * result + (y == null ? 0 : y.hashCode());
			result = prime * result + (z == null ? 0 : z.hashCode());
			return result;
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
		
		public Atom toAtom() {
			AtomImpl a = new AtomImpl();
			a.setCoords(new double[] {x,y,z});
			return a;
		}
	}

	private static final long serialVersionUID = -1523140268972093071L;

	private SymmetryAxisVector axis;
	private SymmetryAxisVector origin;
	private Float orthogonal;
	private Float parallel;
	private Float angle;

	public CensusAxis() {
	}

	public CensusAxis(Float orthogonal, Float screw, Float angle, SymmetryAxisVector axis, SymmetryAxisVector offset) {
		super();
		this.orthogonal = orthogonal;
		this.parallel = screw;
		this.angle = angle;
		this.axis = axis;
		this.origin = offset;
	}

	public CensusAxis(RotationAxis rot) {
		angle = (float) rot.getAngle();
		axis = new SymmetryAxisVector();
		origin = new SymmetryAxisVector();
		if (rot.getRotationAxis() != null) {
			axis.x = (float) rot.getRotationAxis().getX();
			axis.y = (float) rot.getRotationAxis().getY();
			axis.z = (float) rot.getRotationAxis().getZ();
		}
		if (rot.getRotationPos() != null) {
			origin.x = (float) rot.getRotationPos().getX();
			origin.y = (float) rot.getRotationPos().getY();
			origin.z = (float) rot.getRotationPos().getZ();
		}
		if (rot.getScrewTranslation() != null && rot.getRotationAxis() != null) {
			parallel = (float) (Calc.amount(rot.getScrewTranslation()) / Calc.amount(rot.getRotationAxis()));
		}
		if (rot.getOtherTranslation() != null && rot.getRotationAxis() != null) {
			orthogonal = (float) (Calc.amount(rot.getOtherTranslation()) / Calc.amount(rot.getRotationAxis()));
		}
	}
	
	public RotationAxis toRotationAxis() {
		Atom axisAtom = axis.toAtom();
		Atom originAtom = origin.toAtom();
		RotationAxis rotAxis = new RotationAxis(axisAtom,originAtom, angle);
		return rotAxis;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		CensusAxis other = (CensusAxis) obj;
		if (axis == null) {
			if (other.axis != null) return false;
		} else if (!axis.equals(other.axis)) return false;
		if (origin == null) {
			if (other.origin != null) return false;
		} else if (!origin.equals(other.origin)) return false;
		if (orthogonal == null) {
			if (other.orthogonal != null) return false;
		} else if (!orthogonal.equals(other.orthogonal)) return false;
		if (parallel == null) {
			if (other.parallel != null) return false;
		} else if (!parallel.equals(other.parallel)) return false;
		if (angle == null) {
			if (other.angle != null) return false;
		} else if (!angle.equals(other.angle)) return false;
		return true;
	}

	public Double evaluateEpsilon(int order) {
		if (order < 2) return null;
		final double gamma = 2 * Math.PI / order;
		double theta = Math.abs(angle);
		double minEpsilon = Math.abs(theta - gamma);
		double delta = gamma;
		for (int k = 1; delta <= theta; k++, delta += gamma) {
			if (order % k != 0) continue;
			double epsilon = Math.abs(theta - delta);
			if (epsilon < minEpsilon) minEpsilon = epsilon;
		}
		return minEpsilon;
	}

	public SymmetryAxisVector getAxis() {
		return axis;
	}

	public SymmetryAxisVector getOrigin() {
		return origin;
	}

	public Float getOrthogonal() {
		return orthogonal;
	}

	public Float getParallel() {
		return parallel;
	}

	public Float getAngle() {
		return angle;
	}

	public int guessOrder() {
		return guessOrder(1.0 * Math.PI / 180, 8);
	}

	/**
	 * Guesses an order of rotational symmetry from {@link #getTheta() theta}.
	 * 
	 * @param threshold
	 *            The maximal deviation of theta from the expected angle
	 *            (2π/order) required
	 * @param maxOrder
	 *            The maximum order of symmetry to consider
	 * @return The order between 2 and {@code maxOrder} that results in the
	 *         lowest deviation |θ − 2π/order|, or {@code 1} if no deviation
	 *         less than {@code threshold} is found
	 */
	public int guessOrder(double threshold, int maxOrder) {
		return guessOrder(angle, threshold, maxOrder);
	}
	
	public static int guessOrder(double angle, double threshold, int maxOrder) {
		double bestDelta = threshold;
		int bestOrder = 1;
		for (int order = 2; order < maxOrder; order++) {
			double delta = Math.abs(2 * Math.PI / order - angle);
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
		result = prime * result + (axis == null ? 0 : axis.hashCode());
		result = prime * result + (origin == null ? 0 : origin.hashCode());
		result = prime * result + (orthogonal == null ? 0 : orthogonal.hashCode());
		result = prime * result + (parallel == null ? 0 : parallel.hashCode());
		result = prime * result + (angle == null ? 0 : angle.hashCode());
		return result;
	}

	public void setAxis(SymmetryAxisVector axis) {
		this.axis = axis;
	}

	public void setOrigin(SymmetryAxisVector origin) {
		this.origin = origin;
	}

	public void setOrthogonal(Float orthogonal) {
		this.orthogonal = orthogonal;
	}

	public void setParallel(Float parallel) {
		this.parallel = parallel;
	}

	public void setAngle(Float angle) {
		this.angle = angle;
	}

	@Override
	public String toString() {
		return "SymmetryAxis [orthogonal=" + orthogonal + ", screw=" + parallel + ", theta=" + angle + ", axis=" + axis
				+ ", offset=" + origin + "]";
	}

}
