package org.biojava3.structure.align.symm.quaternary;

import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

public class Octahedron implements Polyhedron {
	private static int[] lineLoop1 = {2,4,3,5,2,1,3,0,5,1,4,0,2};
	private double cirumscribedRadius = 1.0;

	/**
	 * Returns the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @return the cirumscribedRadius
	 */
	public double getCirumscribedRadius() {
		return cirumscribedRadius;
	}

	/**
	 * Set the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @param cirumscribedRadius the cirumscribedRadius to set
	 */
	public void setCirumscribedRadius(double cirumscribedRadius) {
		this.cirumscribedRadius = cirumscribedRadius;
	}
	/**
	 * Returns the radius of an inscribed sphere, that is tangent to each 
	 * of the octahedron's faces
	 * @return the inscribedRadius
	 */
	public double getInscribedRadius() {
		double side = getSideLengthFromCircumscribedRadius(cirumscribedRadius);
		return getInscribedRadiusFromSideLength(side);
	}

	/**
	 * Sets the radius of an inscribed sphere, that is tangent to each 
	 * of the octahedron's faces
	 * @param inscribedRadius the inscribedRadius to set
	 */
	public void setInscribedRadius(double radius) {
		double side = getSideLengthFromInscribedRadius(radius);
		this.cirumscribedRadius = getCircumscribedRadiusFromSideLength(side);
	}

	/**
	 * Returns the radius of a sphere, that is tangent to each 
	 * of the octahedron's edges
	 *
	 * @return the midRadius
	 */
	public double getMidRadius() {
		double side = getSideLengthFromCircumscribedRadius(cirumscribedRadius);
		return getMiddleRadiusFromSideLength(side);
	}

	/**
	 * Sets the radius of radius of a sphere, that is tangent to each 
	 * of the octahedron's edges
	 * @param midRadius the midRadius to set
	 */
	public void setMidRadius(double radius) {
		double side = getSideLengthFromMiddleRadius(radius);
		this.cirumscribedRadius = getCircumscribedRadiusFromSideLength(side);
	}

	/**
	 * Returns the vertices of an n-fold polygon of given radius and center	
	 * @param n
	 * @param radius
	 * @param center
	 * @return
	 */ 
	public Point3d[] getVertices() {
		Point3d[] octahedron = new Point3d[6];
	    octahedron[0] = new Point3d(-cirumscribedRadius, 0, 0);
	    octahedron[1] = new Point3d( cirumscribedRadius, 0, 0);
	    octahedron[2] = new Point3d(0, -cirumscribedRadius, 0);
	    octahedron[3] = new Point3d(0,  cirumscribedRadius, 0);
	    octahedron[4] = new Point3d(0, 0, -cirumscribedRadius);
	    octahedron[5] = new Point3d(0, 0,  cirumscribedRadius);

		return octahedron;
	};
	
	public List<int[]> getLineLoops() {
		return Arrays.asList(lineLoop1);
	}
	
	private static double getSideLengthFromInscribedRadius(double radius) {
		return radius * 6 / Math.sqrt(6);
	}
	
	private static double getInscribedRadiusFromSideLength(double sideLength) {
		return sideLength / 6 * Math.sqrt(6);
	}
	
	private static double getSideLengthFromMiddleRadius(double radius) {
		return radius * 2;
	}
	
	private static double getMiddleRadiusFromSideLength(double sideLength) {
		return sideLength / 2;
	}
	
	private static double getSideLengthFromCircumscribedRadius(double radius) {
		return radius * 2 / Math.sqrt(2);
	}
	
	private static double getCircumscribedRadiusFromSideLength(double sideLength) {
		return sideLength / 2 * Math.sqrt(2);
	}
}
