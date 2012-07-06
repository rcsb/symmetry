package org.biojava3.structure.align.symm.quaternary;

import java.util.Random;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.GMatrix;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Structure;

public class AxisTransformation {
	private static Vector3d X_AXIS = new Vector3d(1,0,0);
	private static Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static Vector3d Z_AXIS = new Vector3d(0,0,1);
	
	public static Matrix4d getTransformation(Structure structure, Subunits subunits, RotationGroup rotationGroup) {
		if (subunits.getSubunitCount() == 0) {
			Matrix4d m = new Matrix4d();
			m.setIdentity();
			return m;
		}
		if (rotationGroup.getPointGroup().equals("C1")) {
			return getTransformationByInertiaAxes(subunits);
		} else {
	     	return getTransformationBySymmetryAxes(subunits, rotationGroup);
		}
	}
	
	public static String getJmolTransformation(Matrix4d matrix) {
		Quat4d q = new Quat4d();
		matrix.get(q);
		return "rotate quaternion {" + jMolFloat(q.x) + " " + jMolFloat(q.y) + " " + jMolFloat(q.z) + " " + jMolFloat(q.w) + "}";
	}

	private static Matrix4d getTransformationBySymmetryAxes(Subunits subunits, RotationGroup rotationGroup) {
		
		// orientation of subunit symmetry axes at the centroid of the subunits
		Point3d[] refPoints = new Point3d[3];
		
		// first reference point is centroid of subunits
		refPoints[0] = new Point3d(subunits.getCentroid());
		
		// second reference point is along the principal axis
		Vector3d principalAxis = getPrincipalAxis(rotationGroup);		
		refPoints[1] = new Point3d(subunits.getCentroid());
		refPoints[1].add(principalAxis);
		
		// third reference point is along an orthogonal rotation axis (if available), otherwise, an orthogonal vector
		// to a subunit center is used.
		
		// if rotation group has orthogonal axis, use it for alignment
		Vector3d orthogonalAxis = getOrthogonalAxis(rotationGroup);
		if (orthogonalAxis == null) {
			System.out.println("Ortho is null");
			// need a consistent reference point in case subunits are scrambled in multimers (i.e. use largest subunit)?
			orthogonalAxis = new Vector3d();
			orthogonalAxis.sub(subunits.getOriginalCenters().get(0), refPoints[0]);
		}
		orthogonalAxis.normalize();
		orthogonalAxis.cross(principalAxis, orthogonalAxis);
		orthogonalAxis.normalize();
		
		refPoints[2] = new Point3d(subunits.getCentroid());
		refPoints[2].add(orthogonalAxis);
		
		// check if subunits are already co-linear with principal axis, for example 1A6D, bioassembly 1
		double dp =  Z_AXIS.dot(principalAxis);
		if (Math.abs(dp) > 0.99999) {
			System.out.println("Axis vs. z: " + dp);
			System.out.println("Angle with Y axis: " + Math.toDegrees(X_AXIS.angle(orthogonalAxis)));
			double angle = X_AXIS.angle(orthogonalAxis);
			AxisAngle4d aa = new AxisAngle4d(Z_AXIS, angle);
			Matrix4d m = new Matrix4d();
			m.set(aa);
			return m;
		}
		
		//  y,z axis centered at the centroid of the subunits
		Point3d[] coordPoints = new Point3d[3];
		coordPoints[0] = new Point3d(subunits.getCentroid());
		coordPoints[1] = new Point3d(subunits.getCentroid());
		coordPoints[1].add(Z_AXIS);
		coordPoints[2] = new Point3d(subunits.getCentroid());
		coordPoints[2].add(X_AXIS);
		
		// align principal axis with z axis and perpendicular axis with y axis
		Matrix4d matrix = SuperPosition.superposeWithTranslation(refPoints, coordPoints);

		return matrix;
	}

	private static Vector3d getPrincipalAxis(RotationGroup rotationGroup) {
		Rotation rotation = rotationGroup.getRotation(0);
		AxisAngle4d axisAngle = rotation.getAxisAngle();
		Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		v.normalize();
		return v;
	}
	
	private static Vector3d getOrthogonalAxis(RotationGroup rotationGroup) {
		// find first axis that is orthogonal to principal axis (direction = 1)
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
			if (rotationGroup.getRotation(i).getDirection() == 1) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				v.normalize();
				return v;
			}
		}
		return null;
	}
	
	private static Matrix4d getTransformationByInertiaAxes(Subunits subunits) {
		System.out.println("getTransformationByInertiaAxes");
		Vector3d[] inertiaVectors = subunits.getMomentsOfInertia().getPrincipalAxes();
		
		// orientation of subunit inertia vectors center at the centroid of the subunits
		Point3d[] refPoints = new Point3d[inertiaVectors.length];
		for (int i = 0; i < inertiaVectors.length; i++) {
			refPoints[i] = new Point3d(inertiaVectors[i]);
			refPoints[i].add(subunits.getCentroid());
		}
	
		// x,y,z axis center at the centroid of the subunits
		Point3d[] coordPoints = new Point3d[3];
		coordPoints[0] = new Point3d(X_AXIS);
		coordPoints[0].add(subunits.getCentroid());
		coordPoints[1] = new Point3d(Y_AXIS);
		coordPoints[1].add(subunits.getCentroid());
		coordPoints[2] = new Point3d(Z_AXIS);
		coordPoints[2].add(subunits.getCentroid());

		// align inertia axis with x,y,z axis
		System.out.println("Superposition");
		Matrix4d matrix = SuperPosition.superposeWithTranslation(refPoints, coordPoints);
         System.out.println("return");
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
	
	
	public static String getSymmetryAxesJmol(Subunits subunits, RotationGroup rotationGroup) {
		double radius = subunits.getMomentsOfInertia().getRadiusOfGyration() * 1.2;
		StringBuilder s = new StringBuilder();
		int n = rotationGroup.getOrder();
		float diameter = 1;
		String color = "white";

		for (int i = 0; i < n; i++) {
			Rotation rotation = rotationGroup.getRotation(i);
			
			// don't draw redundant n-fold rotations around principal axis
			if (i > 0 && rotation.getDirection() == 0) {
				continue;
			}
			
			AxisAngle4d axisAngle = rotation.getAxisAngle();
			Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
			v.normalize();
		
			Point3d p1 = new Point3d(v);
			p1.scaleAdd(-radius, subunits.getCentroid());

			Point3d p2 = new Point3d(v);
			p2.scaleAdd(radius, subunits.getCentroid());

			s.append("draw");
			s.append(" c");
			s.append(i);
			s.append(" cylinder {");
			s.append(jMolFloat(p1.x));
			s.append(" ");
			s.append(jMolFloat(p1.y));
			s.append(" ");
			s.append(jMolFloat(p1.z));
			s.append("} {");
			s.append(jMolFloat(p2.x));
			s.append(" ");
			s.append(jMolFloat(p2.y));
			s.append(" ");
			s.append(jMolFloat(p2.z));
			s.append("} diameter ");
			s.append(diameter);
			s.append(" color ");
			s.append(color);
			s.append(";");

			diameter = 0.1f;
			color = "gray";
			s.append(getPolygonJmol(i, p1, p2, v, rotation.getFold()));
			s.append(getPolygonJmol(i + n, p2, p1, v, rotation.getFold()));
		}
	
		return s.toString();
	}
	
	private static String getPolygonJmol(int index, Point3d p1, Point3d p2, Vector3d axis, int n) {
		StringBuilder s = new StringBuilder();
		s.append("draw p");
		s.append(index);
		s.append(" ");
		s.append("polygon");
		s.append(" ");
		s.append(n+1);
		s.append(" ");
		
		s.append("{");
    	s.append(jMolFloat(p1.x));
    	s.append(" ");
    	s.append(jMolFloat(p1.y));
    	s.append(" ");
    	s.append(jMolFloat(p1.z));
      	s.append("}");

        Vector3d[] v = getPolygon(axis, p1, n);
        // create vertex list
        for (int i = 0; i < n; i++) {
        	s.append("{");
        	s.append(jMolFloat(v[i].x));
        	s.append(" ");
        	s.append(jMolFloat(v[i].y));
        	s.append(" ");
        	s.append(jMolFloat(v[i].z));
          	s.append("}");
        }
        
        // create face list
        s.append(" ");
        s.append(n);
        s.append(" ");
        for (int i = 1; i <= n; i++) {
        	s.append("[");
        	s.append(0);
        	s.append(" ");
        	s.append(i);
        	s.append(" ");
        	if (i < n) {
              s.append(i+1);
        	} else {
              s.append(1);
        	}
        	s.append(" ");
        	s.append(7);
           	s.append("]");
        }
        s.append(" mesh");
		s.append(" color red;");
   
		return s.toString();
	}
	
	private static Vector3d[] getPolygon(Vector3d axis, Point3d point, int n) {
		Vector3d perp = getPerpendicularVector(axis);

		AxisAngle4d axisAngle = new AxisAngle4d(axis, 0);
		Vector3d[] vectors = new Vector3d[n];
		System.out.println("Drawing polygon: " + n);
		for (int i = 0; i < n; i++) {
			axisAngle.setAngle(i * 2 * Math.PI/n);
			vectors[i] = new Vector3d(perp);	
			Matrix4d m = new Matrix4d();
			m.set(axisAngle);
			m.transform(vectors[i]);
			vectors[i].add(point);
		}
		return vectors;
	}

	private static Vector3d getPerpendicularVector(Vector3d axis) {
		Random r = new Random();
		Vector3d perp = new Vector3d(r.nextDouble(), r.nextDouble(), r.nextDouble());
		perp.cross(perp, axis);
		perp.normalize();
		return perp;
	}
	
	private static float jMolFloat(double f) {
		if (Math.abs(f) < 1.0E-5) {
			return 0.0f;
		}
		return (float)f;
	}
}
