package org.biojava3.structure.align.symm.quaternary;

import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

public class RectangularPrism implements Polyhedron {
	private static int[] lineLoop1 = {0,1,2,3,0,4,5,6,7,4};
	private static int[] lineLoop2 = {1,5};
	private static int[] lineLoop3 = {2,6};
	private static int[] lineLoop4 = {3,7};
	private double length = 1.0;
	private double width = 1.0;
	private double height = 1.0;

	public RectangularPrism(double length, double width, double height) {
		System.out.println("rectangle: " + length + " " + width + " " + height);
		this.length = length;
		this.width = width;
		this.height = height;
	}
	
	/**
	 * Returns the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @return the cirumscribedRadius
	 */
	public double getLength() {
		return length;
	}

	/**
	 * Returns the radius of an inscribed sphere, that is tangent to each 
	 * of the octahedron's faces
	 * @return the inscribedRadius
	 */
	public double getWidth() {
		return width;
	}

	/**
	 * Returns the radius of a sphere, that is tangent to each 
	 * of the octahedron's edges
	 *
	 * @return the midRadius
	 */
	public double getHeight() {
        return height;
	}

	/**
	 * Returns the vertices of an n-fold polygon of given radius and center	
	 * @param n
	 * @param radius
	 * @param center
	 * @return
	 */ 
	public Point3d[] getVertices() {
		Point3d[] vertices = new Point3d[8];
	    vertices[0] = new Point3d(-width/2, -height/2, length/2);
	    vertices[1] = new Point3d(-width/2, height/2, length/2);
	    vertices[2] = new Point3d(width/2, height/2, length/2);
	    vertices[3] = new Point3d(width/2, -height/2, length/2);
	    vertices[4] = new Point3d(-width/2, -height/2, -length/2);
	    vertices[5] = new Point3d(-width/2, height/2, -length/2);
	    vertices[6] = new Point3d(width/2, height/2, -length/2);
	    vertices[7] = new Point3d(width/2, -height/2, -length/2);

		return vertices;
	};
	
	public List<int[]> getLineLoops() {
		return Arrays.asList(lineLoop1, lineLoop2, lineLoop3, lineLoop4);
	}
	
}
