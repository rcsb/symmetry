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
public final class JmolGeometricObjects {
	
	public static String getJmolWirePolyhedron(int index, Polyhedron p, Matrix4d transformation) {
		String color = "darkorange";
		StringBuilder s = new StringBuilder();

		Point3d[] vertices = p.getVertices();
		for (int i = 0; i < vertices.length; i++) {
			transformation.transform(vertices[i]);
		}

		for (int[] lineLoop: p.getLineLoops()) {
			s.append("draw l");
			s.append(index++);
			s.append(" line ");
			for (int i: lineLoop) {
				s.append(getJmolPoint(vertices[i]));
			}
			s.append(" width 2.0");
			s.append(" color ");
			s.append(color);
			s.append(";");
		}

		return s.toString();
	}
	
	private static String getJmolPoint(Point3d point) {
		StringBuilder s = new StringBuilder();
		s.append("{");
		s.append(jMolFloat(point.x));
		s.append(" ");
		s.append(jMolFloat(point.y));
		s.append(" ");
		s.append(jMolFloat(point.z));
		s.append("}");
		return s.toString();
	}

	public static String getSymmetryAxis(int index, int direction, String pointGroup, int n, Vector3d principalAxis, Vector3d referenceAxis, double radius, float diameter, String color, Point3d center, Vector3d axis) {
		if ((pointGroup.startsWith("C") ||pointGroup.startsWith("D")) && direction == 1) {
			Prism p = new Prism(n);
			p.setInscribedRadius(radius);
			radius = p.getCirumscribedRadius();
	    }
		if (pointGroup.equals("T")) {
			Tetrahedron t = new Tetrahedron();
			t.setMidRadius(radius);
			radius = t.getCirumscribedRadius();
		}
		if (pointGroup.equals("O")) {
			Octahedron o = new Octahedron();
			o.setMidRadius(radius);
			radius = o.getCirumscribedRadius();
		}
		
		Point3d p1 = new Point3d(axis);
		p1.scaleAdd(-radius, center);

		Point3d p2 = new Point3d(axis);
		p2.scaleAdd(radius, center);

		StringBuilder s = new StringBuilder();
		s.append("draw");
		s.append(" c");
		s.append(index);
		s.append(" cylinder ");
		s.append(getJmolPoint(p1));
		s.append(getJmolPoint(p2));
		s.append(" diameter ");
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
			s.append(getLineJmol(index + n + 1, p2, axis, principalAxis, referenceAxis, color));
		} else if (n > 2) {
			double polygonRadius = 3;
			s.append(getPolygonJmol(index, p1, principalAxis, referenceAxis, axis, n, color, polygonRadius));
			s.append(getPolygonJmol(index + n + 1, p2, principalAxis, referenceAxis, axis, n, color, polygonRadius));
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

		Vector3d[] vertexes = getPolygonVertices(principalAxis, axis, referenceAxis, center, 2, 3);
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

		s.append(" width 1.0 ");
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
