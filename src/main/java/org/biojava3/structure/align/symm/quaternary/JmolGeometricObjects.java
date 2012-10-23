/**
 * 
 */
package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * @author Peter
 *
 */
public final class JmolGeometricObjects {

	public static String getJmolPlainPolygon(int index, int n, double height, double radius, Matrix4d transformation) {
		String color = "orange";
		StringBuilder s = new StringBuilder();

		System.out.println("PlainPolygon: " + n + " " + height + " " + radius);
		System.out.println(transformation);
		Point3d center1 = new Point3d(0, 0, -height/2);
		Point3d[] polygon1 = getPolygonVertices(n, radius, center1);
		for (int i = 0; i < polygon1.length; i++) {
			transformation.transform(polygon1[i]);
		}

		s.append("draw l");
		s.append(index);
		s.append(" ");
		s.append("line"); 	
		s.append(" ");

		// create vertex list
		for (Point3d v: polygon1) {
			s.append("{");
			s.append(jMolFloat(v.x));
			s.append(" ");
			s.append(jMolFloat(v.y));
			s.append(" ");
			s.append(jMolFloat(v.z));
			s.append("}");
		}

		s.append("{");
		s.append(jMolFloat(polygon1[0].x));
		s.append(" ");
		s.append(jMolFloat(polygon1[0].y));
		s.append(" ");
		s.append(jMolFloat(polygon1[0].z));
		s.append("}");
		s.append(" color ");
		s.append(" ");
		s.append(color);
		s.append(";");

		Point3d center2 = new Point3d(0,0,height/2);
		System.out.println("getJmolPlainPolygon radius: " + radius);
		Point3d[] polygon2 = getPolygonVertices(n, radius, center2);
		for (int i = 0; i < polygon2.length; i++) {
			transformation.transform(polygon2[i]);
		}

		s.append("draw l");
		s.append(index+1);
		s.append(" ");
		s.append("line"); 	
		s.append(" ");

		for (Point3d v: polygon2) {
			s.append("{");
			s.append(jMolFloat(v.x));
			s.append(" ");
			s.append(jMolFloat(v.y));
			s.append(" ");
			s.append(jMolFloat(v.z));
			s.append("}");
		}

		s.append("{");
		s.append(jMolFloat(polygon2[0].x));
		s.append(" ");
		s.append(jMolFloat(polygon2[0].y));
		s.append(" ");
		s.append(jMolFloat(polygon2[0].z));
		s.append("}");
		s.append(" color ");
		s.append(" ");
		s.append(color);
		s.append(";");

		for (int i = 0; i < polygon1.length; i++) {
			Point3d v1 = polygon1[i];
			Point3d v2 = polygon2[i];

			s.append("draw l");
			s.append(index+2+i);
			s.append(" line ");

			s.append("{");
			s.append(jMolFloat(v1.x));
			s.append(" ");
			s.append(jMolFloat(v1.y));
			s.append(" ");
			s.append(jMolFloat(v1.z));
			s.append("}");
			s.append("{");
			s.append(jMolFloat(v2.x));
			s.append(" ");
			s.append(jMolFloat(v2.y));
			s.append(" ");
			s.append(jMolFloat(v2.z));
			s.append("}");
			s.append(" color ");
			s.append(color);
			s.append(";");

		}
		return s.toString();
	}

	public static String getSymmetryAxis(int index, int n, Vector3d principalAxis, Vector3d referenceAxis, double radius, float diameter, String color, Point3d center, Vector3d axis) {
		Point3d p1 = new Point3d(axis);
		p1.scaleAdd(-radius, center);

		Point3d p2 = new Point3d(axis);
		p2.scaleAdd(radius, center);

		StringBuilder s = new StringBuilder();
		s.append("draw");
		s.append(" c");
		s.append(index);
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
		//		s.append(" color translucent ");
		s.append(" color ");
		s.append(color);
		s.append(";");

		p1 = new Point3d(axis);
		p1.scaleAdd(-radius*0.95, center);

		p2 = new Point3d(axis);
		p2.scaleAdd(radius*0.95, center);

		if (n == 2) {
			s.append(getLineJmol(index, p1, axis, principalAxis, referenceAxis, color));
			s.append(getLineJmol(index + 50, p2, axis, principalAxis, referenceAxis, color));
		} else {
			double polygonRadius = 2;
			s.append(getPolygonJmol(index, p1, principalAxis, referenceAxis, axis, n, color, polygonRadius));
			s.append(getPolygonJmol(index + n, p2, principalAxis, referenceAxis, axis, n, color, polygonRadius));
		}

		return s.toString();
	}
	
	private static String getLineJmol(int index, Point3d center, Vector3d axis, Vector3d principalAxis, Vector3d referenceAxis, String color) {
		StringBuilder s = new StringBuilder();
		s.append("draw l");
		s.append(index);
		s.append(" ");
		s.append("line");
		s.append(" ");

		Vector3d[] vertexes = getPolygonVertices(principalAxis, axis, referenceAxis, center, 2, 2);
		// create vertex list
		for (Vector3d v: vertexes) {
			s.append("{");
			s.append(jMolFloat(v.x));
			s.append(" ");
			s.append(jMolFloat(v.y));
			s.append(" ");
			s.append(jMolFloat(v.z));
			s.append("}");
		}

		s.append(" width 0.5 ");
		s.append(" color ");
		s.append(color);
		s.append(";");

		return s.toString();
	}

	private static String getPolygonJmol(int index, Point3d center, Vector3d principalAxis, Vector3d referenceAxis, Vector3d axis, int n, String color, double radius) {
		StringBuilder s = new StringBuilder();
		s.append("draw p");
		s.append(index);
		s.append(" ");
		s.append("polygon");
		s.append(" ");
		s.append(n+1); 
		s.append(" ");

		s.append("{");
		s.append(jMolFloat(center.x));
		s.append(" ");
		s.append(jMolFloat(center.y));
		s.append(" ");
		s.append(jMolFloat(center.z));
		s.append("}");

		Vector3d[] vertexes = getPolygonVertices(principalAxis, axis, referenceAxis, center, n, radius);
		//	System.out.println("AxisTransformation: polygon: " + Arrays.toString(vertexes));
		// create vertex list
		for (Vector3d v: vertexes) {
			s.append("{");
			s.append(jMolFloat(v.x));
			s.append(" ");
			s.append(jMolFloat(v.y));
			s.append(" ");
			s.append(jMolFloat(v.z));
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
		//	s.append(" color translucent ");
		s.append(" color ");
		s.append(color);
		s.append(";");

		return s.toString();
	}
	private static Vector3d[] getPolygonVertices(Vector3d principalAxis, Vector3d axis, Vector3d referenceAxis, Point3d center, int n, double radius) {
		Vector3d perp = new Vector3d(axis);
		//	System.out.println("AxisTransformation: radius: " + radius);
		//	System.out.println("AxisTransformation: center: " + center);
		//	System.out.println("AxisTransformation: principal axis: " + principalAxis);
		//	System.out.println("AxisTransformation: reference axis: " + referenceAxis);
		//	System.out.println("AxisTransformation: axis: " + axis);
		// if axis coincides with principal axis, use the reference axis to orient polygon
		// TODO need to check alignment with y-axis for all cases of symmetry
		if (Math.abs(axis.dot(principalAxis)) > 0.9) {
			perp.set(referenceAxis);
		} else {
			perp.cross(perp, principalAxis);
		}
		//	System.out.println("AxisTransformation: perp. axis: " + referenceAxis);
		perp.scale(radius);		

		AxisAngle4d axisAngle = new AxisAngle4d(axis, 0);
		Vector3d[] vectors = new Vector3d[n];
		Matrix4d m = new Matrix4d();

		for (int i = 0; i < n; i++) {
			axisAngle.angle = i * 2 * Math.PI/n;
			vectors[i] = new Vector3d(perp);		
			m.set(axisAngle);
			m.transform(vectors[i]);
			vectors[i].add(center);
		}
		return vectors;
	}
	/**
	 * Returns the vertices of an n-fold polygon of given radius and center	
	 * @param n
	 * @param radius
	 * @param center
	 * @return
	 */
	private static Point3d[] getPolygonVertices(int n, double radius, Point3d center) {
		Point3d[] vertices = new Point3d[n];
		Matrix3d m = new Matrix3d();

		for (int i = 0; i < n; i++) {
			vertices[i] = new Point3d(0, radius, 0);			
			m.rotZ(i*2*Math.PI/n);
			m.transform(vertices[i]);
			vertices[i].add(center);
		}
		return vertices;
	}
	/**
	 * Returns a lower precision floating point number for Jmol
	 * @param f
	 * @return
	 */
	private static float jMolFloat(double f) {
		if (Math.abs(f) < 1.0E-7) {
			return 0.0f;
		}
		return (float)f;
	}
}
