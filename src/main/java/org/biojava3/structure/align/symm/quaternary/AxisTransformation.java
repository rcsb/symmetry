package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.GMatrix;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Structure;

public class AxisTransformation {
	private static double scaleFactor = 1;
	
	public static Matrix4d getTransformation(Structure structure, Subunits subunits, RotationGroup rotationGroup) {
		if (rotationGroup.getPointGroup().equals("C1")) {
			Matrix4d identity = new Matrix4d();
			identity.setIdentity();
			return identity;
		}
		Point3d[] refPoints = new Point3d[3];
		// first reference point is centroid of subunits
		refPoints[0] = new Point3d(subunits.getCentroid());
//		System.out.println("r0: " + refPoints[0]);
		
		// second reference point is along the principal axis
		Rotation rotation = rotationGroup.getRotation(0);
		AxisAngle4d axisAngle = rotation.getAxisAngle();
		Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		v.normalize();
		refPoints[1] = new Point3d(v);
		refPoints[1].scaleAdd(scaleFactor, refPoints[0]);
//		System.out.println("r1: " + refPoints[1]);

		
		// third reference point is along the line from the centroid of the structure to the centroid of the first subunit
		Vector3d ortho = new Vector3d();
		ortho.sub(subunits.getOriginalCenters().get(0), refPoints[0]);
		ortho.normalize();
		ortho.cross(v, ortho);
		ortho.normalize();
		
		refPoints[2] = new Point3d(ortho);
		refPoints[2].scaleAdd(scaleFactor, refPoints[0]);
		
		double dp =  v.dot(new Vector3d(0,0,1));
		if (dp < -0.99 || dp > 0.99) {
			System.out.println("Axis vs. z: " + v.dot(new Vector3d(0,0,1)));
			Vector3d yp = new Vector3d(0, 1, 0);
			System.out.println("Angle with Y axis: " + Math.toDegrees(yp.angle(ortho)));
			double angle = yp.angle(ortho);
			AxisAngle4d aa = new AxisAngle4d(new Vector3d(0,0,1), angle);
			Matrix4d m = new Matrix4d();
			m.set(aa);
			return m;
		}
//		System.out.println("r2: " + refPoints[2]);
//		System.out.println("Dot product: " + v.dot(ortho));
		
		// center - y-axis - z-axis are the new coordinates
		Point3d[] coordPoints = new Point3d[3];
		coordPoints[0] = new Point3d(refPoints[0]);
//		System.out.println("c0: " + coordPoints[0]);
		coordPoints[1] = new Point3d(0, 0, 1);
		coordPoints[1].scaleAdd(scaleFactor, refPoints[0]);
//		System.out.println("c1: " + coordPoints[1]);
		coordPoints[2] = new Point3d(0, 1, 0);
		coordPoints[2].scaleAdd(scaleFactor, refPoints[0]);
//		System.out.println("c2: " + coordPoints[2]);
		Matrix4d matrix = SuperPosition.superposeWithTranslation(refPoints, coordPoints);
//		Matrix4d matrix = SuperPosition.superpose(refPoints, coordPoints);
//		System.out.println(matrix);
//		System.out.println("After superposition rmsd: " + SuperPosition.rmsd(refPoints, coordPoints));
//		System.out.println("r0: " + refPoints[0]);
//		System.out.println("r1: " + refPoints[1]);
//		System.out.println("r2: " + refPoints[2]);

		return matrix;
	}
	
	public static double getTrace(Matrix4d matrix) {
		GMatrix m = new GMatrix(4,4);
		m.set(matrix);
		System.out.println("Trace: " + m.trace());
		if (m.trace() <= 0) {
			System.out.println(matrix);
		}
		return m.trace();
	}
	
	public static String getJmolQuat(Matrix4d matrix) {
		Quat4d q = new Quat4d();
		matrix.get(q);
		return "rotate quaternion {" + q.x + " " + q.y + " " + q.z + " " + q.w + "}";
	}
}
