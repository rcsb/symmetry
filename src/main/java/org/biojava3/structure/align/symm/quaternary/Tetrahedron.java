package org.biojava3.structure.align.symm.quaternary;

import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;

public class Tetrahedron implements Polyhedron {
	private static int[] lineLoop1 = {0,1,2,3,0,2};
	private static int[] lineLoop2 = {1,3};
	
	private double circumscribedRadius = 1.0;

	/**
	 * Returns the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @return the cirumscribedRadius
	 */
	public double getCirumscribedRadius() {
		return circumscribedRadius;
	}

	/**
	 * Set the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @param cirumscribedRadius the cirumscribedRadius to set
	 */
	public void setCirumscribedRadius(double cirumscribedRadius) {
		this.circumscribedRadius = cirumscribedRadius;
	}
	/**
	 * Returns the radius of an inscribed sphere, that is tangent to each 
	 * of the tetrahedrons's faces
	 * @return the inscribedRadius
	 */
	public double getInscribedRadius() {
		double side = getSideLengthFromCircumscribedRadius(circumscribedRadius);
		return getInscribedRadiusFromSideLength(side);
	}

	/**
	 * Sets the radius of an inscribed sphere, that is tangent to each 
	 * of the tetrahedron's faces
	 * @param inscribedRadius the inscribedRadius to set
	 */
	public void setInscribedRadius(double radius) {
		double side = getSideLengthFromInscribedRadius(radius);
		this.circumscribedRadius = getCircumscribedRadiusFromSideLength(side);
	}

	/**
	 * Returns the radius of a sphere, that is tangent to each 
	 * of the tetrahedron's edges
	 *
	 * @return the midRadius
	 */
	public double getMidRadius() {
		double side = getSideLengthFromCircumscribedRadius(circumscribedRadius);
		return getMiddleRadiusFromSideLength(side);
	}

	/**
	 * Sets the radius of radius of a sphere, that is tangent to each 
	 * of the tetrahedron's edges
	 * @param midRadius the midRadius to set
	 */
	public void setMidRadius(double radius) {
		double side = getSideLengthFromMiddleRadius(radius);
		this.circumscribedRadius = getCircumscribedRadiusFromSideLength(side);
	}

	/**
	 * Returns the vertices of an n-fold polygon of given radius and center	
	 * @param n
	 * @param radius
	 * @param center
	 * @return
	 */ 
	public  Point3d[] getVertices() {
		double a = circumscribedRadius/Math.sqrt(2);
	    Point3d[] tetrahedron = new Point3d[4];
		tetrahedron[0] = new Point3d(-circumscribedRadius, 0, -a);
		tetrahedron[1] = new Point3d( circumscribedRadius, 0, -a);
		tetrahedron[2] = new Point3d(0, -circumscribedRadius, a);
		tetrahedron[3] = new Point3d(0,  circumscribedRadius, a);
		Point3d centroid = SuperPosition.centroid(tetrahedron);
		
		// rotate tetrahedron to align one vertex with the +z axis
		Matrix3d m = new Matrix3d();
		double tetrahedralAngle = Math.acos(-1.0/3.0);
		m.rotX(0.5 * tetrahedralAngle);
		for (Point3d p: tetrahedron) {
			p.sub(centroid);
			m.transform(p);
		}
		return tetrahedron;
	};
	
	public List<int[]> getLineLoops() {
		return Arrays.asList(lineLoop1, lineLoop2);
	}
	
	private static double getSideLengthFromInscribedRadius(double radius) {
		return radius * Math.sqrt(24);
	}
	
	private static double getInscribedRadiusFromSideLength(double sideLength) {
		return sideLength / Math.sqrt(24);
	}
	
	private static double getSideLengthFromMiddleRadius(double radius) {
		return radius * Math.sqrt(8);
	}
	
	private static double getMiddleRadiusFromSideLength(double sideLength) {
		return sideLength / Math.sqrt(8);
	}
	
	private static double getSideLengthFromCircumscribedRadius(double radius) {
		return radius / Math.sqrt(3.0/8.0);
	}
	
	private static double getCircumscribedRadiusFromSideLength(double sideLength) {
		return sideLength * Math.sqrt(3.0/8.0);
	}
	
}
