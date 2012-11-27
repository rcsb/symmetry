/**
 * 
 */
package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * @author Peter
 *
 */
public class SymmetrySymbol {

	public static Point3d[] getPolygon(int n, Vector3d axis, Point3d center, Point3d[] vertices) {
		double radius = 2;
	//	Point3d[] vertices = polyhedron.getVertices();
		double minDist = Double.MAX_VALUE;
		for (Point3d p: vertices) {
			double dSq = p.distanceSquared(center);
			if (dSq > 0.01) {
				minDist = Math.min(p.distanceSquared(center), minDist);
				System.out.println("dSq: " + dSq);
			}
		}
		System.out.println("minDist: " + minDist);

		int index = 0;
		Point3d[] referencePoints = new Point3d[n];
		for (Point3d p: vertices) {
			double dSq = p.distanceSquared(center);
			if (dSq > 0.1 && Math.abs(dSq-minDist) < 0.1) {
				referencePoints[index] = p;
				index++;
				System.out.println("index: " + index);
			}
		}
		
		Point3d[] refPoints = new Point3d[index];
		for (int i = 0; i < refPoints.length; i++) {
			Vector3d v = new Vector3d();
			v.set(referencePoints[i]);
			v.sub(center);
			v.cross(axis, v);
			v.cross(axis, v);
			referencePoints[i].set(v);
			referencePoints[i].scale(radius);
			referencePoints[i].add(center);
		}
		
		return refPoints;
	}
}
