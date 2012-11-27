/**
 * 
 */
package org.biojava3.structure.align.symm.jmolScript;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.align.symm.geometry.Polyhedron;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.Rotation;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;
import org.biojava3.structure.align.symm.quaternary.SymmetrySymbol;

/**
 * @author Peter
 *
 */
public abstract class JmolSymmetryScriptGenerator {
	private static String POLYHEDRON_COLOR = "lawngreen";
	private static String N_FOLD_AXIS_COLOR = "red";
	private static String TWO_FOLD_AXIS_COLOR = "blue";
	private static String THREE_FOLD_AXIS_COLOR = "green";

	private static double AXIS_SCALE_FACTOR = 1.2;
	
	protected AxisTransformation axisTransformation = null;
	protected RotationGroup rotationGroup = null;
	protected Polyhedron polyhedron = null;
	
	public JmolSymmetryScriptGenerator(AxisTransformation axisTransformation, RotationGroup rotationGroup) {
		this.axisTransformation = axisTransformation;
		this.rotationGroup = rotationGroup;
	}

	abstract public int getDefaultZoom();

	public String deleteSymmetryAxes() {return new String();};
	public String deleteInertiaAxes() {return new String();};
	
	public int getOrientationCount() {
		return polyhedron.getViewCount();
	}
	
	public String setOrientation(int index) {	
		StringBuilder s = new StringBuilder();
		s.append(setCentroid());
		
		// calculate  orientation
		Matrix3d m = polyhedron.getViewMatrix(index);
		m.mul(axisTransformation.getRotationMatrix());
		Quat4d q = new Quat4d();
		q.set(m);
		
		// set orientation
		s.append("moveto 4 quaternion {");
		s.append(jMolFloat(q.x));
		s.append(" ");
		s.append(jMolFloat(q.y));
		s.append(" ");
		s.append(jMolFloat(q.z));
		s.append(" ");
		s.append(jMolFloat(q.w));
		s.append("} ");
		s.append(getDefaultZoom());
		s.append(";");
		return s.toString();
	}
	
	public String drawOrientationName(int index, String color) {
		StringBuilder s = new StringBuilder();
		s.append("set echo top center;");
		s.append("color echo ");
		s.append(color);
		s.append(";");
		s.append("font echo 24 sanserif;");
		s.append("echo ");
		s.append(polyhedron.getViewName(index));
		s.append(";");
		return s.toString();
	}

	public String drawPolyhedron() {
		StringBuilder s = new StringBuilder();

		Point3d[] vertices = polyhedron.getVertices();
		Matrix4d reverseTransformation = axisTransformation.calcPrismTransformation();
		for (int i = 0; i < vertices.length; i++) {
			reverseTransformation.transform(vertices[i]);
		}
		
		int index = 100;

		for (int[] lineLoop: polyhedron.getLineLoops()) {
			s.append("draw polyhedron");
			s.append(index++);
			s.append(" line ");
			for (int i: lineLoop) {
				s.append(getJmolPoint(vertices[i]));
			}
			s.append(" width 1.5");
			s.append(" color ");
			s.append(POLYHEDRON_COLOR);
			s.append(" off;");
		}

		return s.toString();
	}
	
	public String hidePolyhedron() {
		return "draw polyhedron* off;";
	}
	
	public String showPolyhedron() {
		return "draw polyhedron* on;";
	}
	
	public String hideAxes() {
		return "draw axes* off;";
	}
	
	public String showAxes() {
		return "draw axes* on;";
	}
	
	public String drawInertiaAxes() {
		StringBuilder s = new StringBuilder();
		Point3d centroid = axisTransformation.calcGeometricCenter();
		Vector3d[] axes = axisTransformation.getPrincipalAxesOfInertia();

		for (int i = 0; i < axes.length; i++) {
			s.append("draw axesInertia");
			s.append(200+i);
			s.append(" ");
			s.append("line");
			s.append(" ");
			Point3d v1 = new Point3d(axes[i]);
			if (i == 0) {
				v1.scale(1.2*axisTransformation.getDimension().y);
			} else if (i == 1) {
				v1.scale(1.2*axisTransformation.getDimension().x);
			} else if (i == 2) {
				v1.scale(1.2*axisTransformation.getDimension().z);
			}
			Point3d v2 = new Point3d(v1);
			v2.negate();
			v1.add(centroid);
			v2.add(centroid);
			s.append(getJmolPoint(v1));
			s.append(getJmolPoint(v2));
			s.append(" width 0.5 ");
			s.append(" color white");
			s.append(" off;");
		}
        return s.toString();
	};
	
	public String drawPointGroupName(String color) {
		StringBuilder s = new StringBuilder();
		s.append("set echo bottom center;");
		s.append("color echo ");
		s.append(color);
		s.append(";");
		s.append("font echo 24 sanserif;");
		s.append("echo Point group ");
		s.append(rotationGroup.getPointGroup());
		s.append(";");
		return s.toString();
	}
	
	public String playOrientations() {
		StringBuilder s = new StringBuilder();
		
		// TODO point group name is not displayed???
		s.append(drawPointGroupName("white"));
		
		// draw polygon
		s.append(drawPolyhedron());
		s.append(showPolyhedron());
		s.append(hidePolyhedron());
	//	s.append(showPolyhedron());
			
		// draw axes
		s.append(drawAxes());
		s.append(showAxes());
		s.append(hideAxes());
//		s.append(showAxes());
		
		// loop over all orientations
		for (int i = 0; i < getOrientationCount(); i++) {
			s.append("echo ;");
			s.append(setOrientation(i));
			s.append(drawOrientationName(i, "white"));
			s.append("delay 4;");
		}
		
		// go back to first orientation
		s.append("echo ;");
		s.append(setOrientation(0));
		s.append(drawOrientationName(0, "white"));
		
		s.append(drawPointGroupName("white"));
		
		return s.toString();
	}

	public String drawAxes() {
		if (rotationGroup.getPointGroup().equals("C1")) {
			return drawInertiaAxes();
		} else {
			return drawSymmetryAxes();
		}
	}
	
	public String drawSymmetryAxes() {
		StringBuilder s = new StringBuilder();

		int n = rotationGroup.getOrder();
		if (n == 0) {
			return s.toString();
		}

		float diameter = 0.5f;
		double radius = 0;
		String color = "";

		List<Rotation> axes = getUniqueAxes();
//		System.out.println("Unique axes: " + axes.size());
		int i = 0;
		for (Rotation r: axes) {
			if (rotationGroup.getPointGroup().startsWith("C") || (rotationGroup.getPointGroup().startsWith("D") && r.getDirection() == 0)) {
				radius =  axisTransformation.getDimension().z; // principal axis uses z-dimension
				color = N_FOLD_AXIS_COLOR;
			} else {
				radius = polyhedron.getCirumscribedRadius();
				// TODO radius of Tetrahedron seems to be too small???
				// System.out.println("Symmetry axis radius: " + radius);
				if (rotationGroup.getPointGroup().equals("T")) {
					radius *= 1.22;
				}
			
				if (r.getFold() == 2) {
					color = TWO_FOLD_AXIS_COLOR;
				} else if (r.getFold() == 3) {
					color = THREE_FOLD_AXIS_COLOR;
				} else {
					color = N_FOLD_AXIS_COLOR;
				}
			}
		

			Point3d center = axisTransformation.calcGeometricCenter();
			AxisAngle4d axisAngle = r.getAxisAngle();
			Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
	//		System.out.println("Unique axes: " + axis + " n: " + r.getFold());
			s.append(getSymmetryAxis(i, i+axes.size(), rotationGroup.getPointGroup(), r.getFold(),  axisTransformation.getPrincipalRotationAxis(), axisTransformation.getRotationReferenceAxis(), radius, diameter, color, center, axis));
	        i++;
		}

		return s.toString();
	}
	
	public String getSymmetryAxis(int i, int j, String pointGroup, int n, Vector3d principalAxis, Vector3d referenceAxis, double radius, float diameter, String color, Point3d center, Vector3d axis) {
		boolean drawPolygon = true;
		
		Point3d p1 = new Point3d(axis);
		p1.scaleAdd(-AXIS_SCALE_FACTOR * radius, center);

		Point3d p2 = new Point3d(axis);
		p2.scaleAdd(AXIS_SCALE_FACTOR * radius, center);

		StringBuilder s = new StringBuilder();
		s.append("draw");
		s.append(" axesSymmetry");
		s.append(i);
		s.append(" cylinder ");
		s.append(getJmolPoint(p1));
		s.append(getJmolPoint(p2));
		s.append(" diameter ");
		s.append(diameter);
		//		s.append(" color translucent ");
		s.append(" color ");
		s.append(color);
		s.append(" off;");

		p1 = new Point3d(axis);
		p1.scaleAdd(-1.01*radius, center);

		p2 = new Point3d(axis);
		p2.scaleAdd(1.01*radius, center);

		if (drawPolygon == true) {
//			if (n == 2) {
//				s.append(getLineJmol(index, p1, axis, principalAxis, referenceAxis, color));
//				s.append(getLineJmol(index + n + 1, p2, axis, principalAxis, referenceAxis, color));
//			} else if (n > 2) {
//				double polygonRadius = 5;
			double polygonRadius = radius * 0.05;
				s.append(getPolygonJmol(i, p1, principalAxis, referenceAxis, axis, n, color, polygonRadius));
				s.append(getPolygonJmol(j, p2, principalAxis, referenceAxis, axis, n, color, polygonRadius));
//			}
		}

		return s.toString();
	}
	
	private static String getPolygonJmol(int index, Point3d center, Vector3d principalAxis, Vector3d referenceAxis, Vector3d axis, int n, String color, double radius) {
		StringBuilder s = new StringBuilder();
		s.append("draw axesSymbol");
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

		if (n == 2) {
	      	s.append(" mesh");
		}
		//	s.append(" color translucent ");
		s.append(" color ");
		s.append(color);
		s.append(" off;");

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
	

    private List<Rotation> getUniqueAxes() {
    	List<Rotation> uniqueRotations = new ArrayList<Rotation>();
    	
    	for (int i = 0, n = rotationGroup.getOrder(); i < n; i++) {
			Rotation rotationI = rotationGroup.getRotation(i);
			AxisAngle4d axisAngleI = rotationI.getAxisAngle();
			Vector3d axisI = new Vector3d(axisAngleI.x, axisAngleI.y, axisAngleI.z);
			
			boolean redundant = false;
			for (Rotation r: uniqueRotations) {
				AxisAngle4d axisAngleJ = r.getAxisAngle();
			    Vector3d axisJ = new Vector3d(axisAngleJ.x, axisAngleJ.y, axisAngleJ.z);
	//		    System.out.println("dot: " + axisI.dot(axisJ));
			    if (Math.abs(axisI.dot(axisJ)) > 0.97) {
			    	redundant = true;
			    	break;
			    }
			}
			
			if (! redundant) {
				uniqueRotations.add(rotationI);
			}
    	}
        return uniqueRotations;
    }

//	public String drawSymmetrySymbolOld(int index) {	
//		String color = "red";
//
////		Matrix3d m = polyhedron.getViewMatrix(index);
//		Matrix3d m1 = new Matrix3d();
//		axisTransformation.getReverseTransformation().getRotationScale(m1);
//		
//		Point3d[] vertices = polyhedron.getVertices();	
//		
//		for (int i = 0; i < vertices.length; i++) {
//			m1.transform(vertices[i]);
//			System.out.println("vi: " + vertices[i]);
//		}
//		
//		double zMax = Double.MIN_VALUE;
//		Point3d center = null;
//		for (Point3d v: vertices) {
//			if (v.z > zMax) {
//				zMax = v.z;
//				center = v;
//			}
//		}
//
//		int n = rotationGroup.getRotation(0).getFold();
//
//		StringBuilder s = new StringBuilder();
//		s.append("draw symmetrySymbol");
//		s.append(" ");
//		s.append("polygon");
//		s.append(" ");
//		s.append(n+1); 
//		s.append(" ");
//
//		s.append("{");
//		s.append(jMolFloat(center.x));
//		s.append(" ");
//		s.append(jMolFloat(center.y));
//		s.append(" ");
//		s.append(jMolFloat(center.z));
//		s.append("}");
//		
//		Vector3d axis = new Vector3d();
//		axis.sub(center, axisTransformation.calcGeometricCenter());
//
//
//		Point3d[] pVertices = SymmetrySymbol.getPolygon(n, axis, center, vertices);
//		// create vertex list
//		for (Point3d v: pVertices) {
//			s.append("{");
//			s.append(jMolFloat(v.x));
//			s.append(" ");
//			s.append(jMolFloat(v.y));
//			s.append(" ");
//			s.append(jMolFloat(v.z));
//			s.append("}");
//		}
//
//		// create face list
//		s.append(" ");
//		s.append(pVertices.length);
//		s.append(" ");
//
//		for (int i = 1; i <= pVertices.length; i++) {
//			s.append("[");
//			s.append(0);
//			s.append(" ");
//			s.append(i);
//			s.append(" ");
//			if (i < n) {
//				s.append(i+1);
//			} else {
//				s.append(1);
//			}
//			s.append(" ");
//			s.append(7);
//			s.append("]");
//		}
//
//		s.append(" mesh");
//		//	s.append(" color translucent ");
//		s.append(" color ");
//		s.append(color);
//		s.append(";");
//
//		return s.toString();
//	}
	
	public static JmolSymmetryScriptGenerator getInstance(AxisTransformation axisTransformation, RotationGroup rotationGroup) {
		String pointGroup = rotationGroup.getPointGroup();
		if (pointGroup.equals("C1")) {
			return new JmolSymmetryScriptGeneratorC1(axisTransformation, rotationGroup);
		} else if (pointGroup.startsWith("C")) {
			return new JmolSymmetryScriptGeneratorCn(axisTransformation, rotationGroup);
		} else if (pointGroup.startsWith("D")) {
			return new JmolSymmetryScriptGeneratorCn(axisTransformation, rotationGroup);
		} else if (pointGroup.equals("T")) {
			return new JmolSymmetryScriptGeneratorT(axisTransformation, rotationGroup);
		} else if (pointGroup.equals("O")) {
			return new JmolSymmetryScriptGeneratorO(axisTransformation, rotationGroup);
		} else if (pointGroup.equals("I")) {
			return new JmolSymmetryScriptGeneratorI(axisTransformation, rotationGroup);
		} 
		
		return null;
	}
	
	private String setCentroid() {
		// calculate center of rotation
		Point3d centroid = axisTransformation.calcGeometricCenter();
			
		// set centroid
		StringBuilder s = new StringBuilder();
		s.append("center {");
		s.append(jMolFloat(centroid.x));
		s.append(" ");
		s.append(jMolFloat(centroid.y));
		s.append(" ");
		s.append(jMolFloat(centroid.z));
		s.append("};");
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
