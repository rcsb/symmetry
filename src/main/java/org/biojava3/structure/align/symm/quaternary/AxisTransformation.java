package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.GMatrix;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

public class AxisTransformation {
	private static Vector3d X_AXIS = new Vector3d(1,0,0);
	private static Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static Vector3d Z_AXIS = new Vector3d(0,0,1);

	private Subunits subunits = null;
	private RotationGroup rotationGroup = null;

	private Matrix4d transformationMatrix = new Matrix4d();
	private Vector3d principalAxis = new Vector3d();
	private Vector3d referenceAxis = new Vector3d();

	public AxisTransformation(Subunits subunits, RotationGroup rotationGroup) {
		this.subunits = subunits;
		this.rotationGroup = rotationGroup;
	}

	public Matrix4d getTransformation() {
		if (subunits == null || subunits.getSubunitCount() == 0) {
			transformationMatrix.setIdentity();
		} else if (rotationGroup.getPointGroup().equals("C1")) {
			transformationMatrix = getTransformationByInertiaAxes();
		} else {
			transformationMatrix = getTransformationBySymmetryAxes();
		}
		return transformationMatrix;
	}

	private Matrix4d getTransformationBySymmetryAxes() {

		// orientation of subunit symmetry axes at the centroid of the subunits
		Point3d[] refPoints = new Point3d[3];

		// first reference point is centroid of subunits
		refPoints[0] = new Point3d(subunits.getCentroid());

		// second reference point is along the principal axis
		principalAxis = getPrincipalRotationAxis();		
		refPoints[1] = new Point3d(subunits.getCentroid());
		refPoints[1].add(principalAxis);

		// third reference point is along an orthogonal rotation axis (if available), otherwise, an orthogonal vector
		// to a subunit center is used.

		// if rotation group has orthogonal axis, use it for alignment
		referenceAxis = getMinorRotationAxis();
		if (referenceAxis == null) {
			referenceAxis = getSubunitReferenceAxis();
		}
		referenceAxis.cross(principalAxis, referenceAxis); // make it perpendicular
		referenceAxis.normalize();

		refPoints[2] = new Point3d(subunits.getCentroid());
		refPoints[2].add(referenceAxis);

		// check if subunits are already co-linear with principal axis, for example 1A6D, bioassembly 1
		double dp =  Z_AXIS.dot(principalAxis);
		if (Math.abs(dp) > 0.99999) {
			//			System.out.println("Axis vs. z: " + dp);
			//			System.out.println("Angle with Y axis: " + Math.toDegrees(Y_AXIS.angle(referenceAxis)));
			double angle = Y_AXIS.angle(referenceAxis);
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
		coordPoints[2].add(Y_AXIS);

		// align principal axis with z axis and perpendicular axis with x axis
		Matrix4d matrix = SuperPosition.superposeWithTranslation(refPoints, coordPoints);

		return matrix;
	}

	/**
	 * Returns vector from largest subunit to centroid of complex
	 * @return
	 */
	private Vector3d getSubunitReferenceAxis() {
		Vector3d orthogonalAxis = new Vector3d();
		int index = subunits.getLargestSubunit();
		orthogonalAxis.sub(subunits.getOriginalCenters().get(index), subunits.getCentroid());
		orthogonalAxis.normalize();
		return orthogonalAxis;
	}

	private double[] getDimensions() {
		double axisScale = 1.2;
		double xMin = Double.MAX_VALUE;
		double xMax = Double.MIN_VALUE;
		double yMin = Double.MAX_VALUE;
		double yMax = Double.MIN_VALUE;
		double zMin = Double.MAX_VALUE;
		double zMax = Double.MIN_VALUE;
		double xyRadiusMaxSq = Double.MIN_VALUE;

		Point3d probe = new Point3d();
		Point3d centroid = subunits.getCentroid();

		for (Point3d[] list: subunits.getTraces()) {
			for (Point3d p: list) {
				probe.set(p);
				probe.sub(centroid);
				transformationMatrix.transform(probe);
				xMin = Math.min(xMin, probe.x);
				xMax = Math.max(xMax, probe.x);
				yMin = Math.min(yMin, probe.y);
				yMax = Math.max(yMax, probe.y);
				zMin = Math.min(zMin, probe.z);
				zMax = Math.max(zMax, probe.z);
				double rSq = probe.x*probe.x + probe.y*probe.y;
				xyRadiusMaxSq = Math.max(xyRadiusMaxSq, rSq);
			}
		}
		double[] dimensions = new double[2];
		dimensions[0] = axisScale * 0.5 * (zMax - zMin); // half of dimension along z-axis (principal rotation axis)
		dimensions[1] = axisScale * 0.5 * Math.max((xMax-xMin), (yMax-yMin)); // max radius for rotation in x-y plane

		return dimensions;
	}

	private Vector3d getPrincipalRotationAxis() {
		Rotation rotation = rotationGroup.getRotation(0); // the rotation around the principal axis is the first rotation
		AxisAngle4d axisAngle = rotation.getAxisAngle();
		Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		v.normalize();
		return v;
	}

	private Vector3d getMinorRotationAxis() {
		// find axis that is not the rotation principal axis (direction = 1)
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

	private int getPrincipalRotationAxisFold() {
		if (rotationGroup.getOrder() == 0) {
			return 1;
		}
		Rotation rotation = rotationGroup.getRotation(0); // the rotation around the principal axis is the first rotation
		return rotation.getFold();
	}

	private int getMinorRotationAxisFold() {
		// find axis that is not the rotation principal axis (direction = 1)
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
			if (rotationGroup.getRotation(i).getDirection() == 1) {
				return rotationGroup.getRotation(i).getFold();
			}
		}
		return 1;
	}

	private Matrix4d getTransformationByInertiaAxes() {
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
		Matrix4d matrix = SuperPosition.superposeWithTranslation(refPoints, coordPoints);
		return matrix;
	}

	public String getJmolTransformation() {
		Quat4d q = new Quat4d();
		transformationMatrix.get(q);
		return "rotate quaternion {" + jMolFloat(q.x) + " " + jMolFloat(q.y) + " " + jMolFloat(q.z) + " " + jMolFloat(q.w) + "};";
	}

	public String getJmolSymmetryAxes() {
		StringBuilder s = new StringBuilder();


		int n = rotationGroup.getOrder();
		if (n == 0) {
			return s.toString();
		}
		double[] dimensions = getDimensions();
		float diameter = 0.5f;
		double radius = 0;
		String color = "red";

		for (int i = 0; i < n; i++) {
			Rotation rotation = rotationGroup.getRotation(i);
			int direction =  rotation.getDirection();

			// don't draw redundant n-fold rotations around principal axis
			if (i > 0 && direction == 0) {
				continue;
			}

			if (direction == 0) {
				radius = dimensions[0]; // principal axis uses z-dimension
				color = "red";
				diameter = 0.5f;
			} else {
				radius = dimensions[1];
				color = "royalblue";
				diameter = 0.25f;
			}

			AxisAngle4d axisAngle = rotation.getAxisAngle();
			Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);

			Point3d p1 = new Point3d(axis);
			p1.scaleAdd(-radius, subunits.getCentroid());

			Point3d p2 = new Point3d(axis);
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
	//		s.append(" color translucent ");
			s.append(" color ");
			s.append(color);
			s.append(";");

			p1 = new Point3d(axis);
			p1.scaleAdd(-radius*0.95, subunits.getCentroid());

			p2 = new Point3d(axis);
			p2.scaleAdd(radius*0.95, subunits.getCentroid());

			if (rotation.getFold() == 2) {
				s.append(getLineJmol(i, p1, axis, color));
				s.append(getLineJmol(i + n, p2, axis, color));
			} else {
				s.append(getPolygonJmol(i, p1, axis, rotation.getFold(), color));
				s.append(getPolygonJmol(i + n, p2, axis, rotation.getFold(), color));
			}
		}

		return s.toString();
	}

	public String getJmolAnimation(int delay) {
		String animation = "";
		String pointGroup = rotationGroup.getPointGroup();
		if (pointGroup.startsWith("C")) {
			animation = getJmolAnimationCyclic(delay);
		} else if (pointGroup.startsWith("D")) {
			animation = getJmolAnimationDihedral(delay);
		} else if (pointGroup.startsWith("O")) {
		    animation = getJmolAnimationOctahedral(delay);
	    } else if (pointGroup.startsWith("I")) {
		    animation = getJmolAnimationIcosahedral(delay);
	    }
		return animation;
	}
	
	private String getJmolAnimationCyclic(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());	
		
		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		int fold  = getPrincipalRotationAxisFold();
		if (fold > 1) {
			// show symmetry axes
			s.append("set echo top center;");
			s.append("echo ");
			s.append(fold);
			s.append("-fold rotation;");
			s.append("color echo red;");
			s.append(getJmolSymmetryAxes());
			s.append("delay ");
			s.append(delay); 
			s.append(";");

			// rotate around principal axis
			fold  = getPrincipalRotationAxisFold();
			float angle = 360.0f/fold;
			s.append("move 0 0 ");
			s.append(angle);
			s.append(" 0 0 0 0 0 3;");
			s.append("delay ");
			s.append(delay);
			s.append(";");
		}
		
		// show Side view 1
		s.append("set echo top center;");
		s.append("echo Side view 1;");
		s.append("color echo white;");
		s.append("move 0 90 0 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		
		// show Back view
		s.append("set echo top center;");
		s.append("echo Back view;");
		s.append("move 0 90 0 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		
		// show Side view 2
		s.append("set echo top center;");
		s.append("echo Side view 2;");
		s.append("move 0 90 0 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		
		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("move 0 90 0 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}
	
	public String getJmolAnimationDihedral(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());

	    // show front view
		s.append("set echo top center;");
		s.append("echo Front view;");
        s.append("color echo white;");
        s.append("font echo 24 sanserif;");
        s.append("set echo bottom center;");
		s.append("echo Point group ");
		s.append(rotationGroup.getPointGroup());
		s.append(";");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		int fold  = getPrincipalRotationAxisFold();
		s.append("set echo top center;");
		s.append("echo ");
		s.append(fold);
		s.append("-fold rotation;");
        s.append("color echo red;");
		s.append(getJmolSymmetryAxes());
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		
		// rotate around principal axis
		float angle = 360.0f/fold;
		s.append("move 0 0 ");
		s.append(angle);
		s.append(" 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// align minor axis with z-axis
		fold  = getMinorRotationAxisFold();
		if (fold > 1) {
			// show side view
			s.append("set echo top center;");
			s.append("echo Side view;");
            s.append("color echo white;");
			// align minor axis with z-axis
			s.append("move 0 90 0 0 0 0 0 0 3;");

			// rotate around minor axis
			fold  = getMinorRotationAxisFold();
			angle = 360.0f/fold;
			s.append("delay ");
			s.append(delay); 
			s.append(";");
			s.append("set echo top center;");
			s.append("echo ");
			s.append(fold);
			s.append("-fold rotation;");
	        s.append("color echo royalblue;");
			s.append("move 0 0 ");
			s.append(angle);
			s.append(" 0 0 0 0 0 3;");
			s.append("delay ");
			s.append(delay);
			s.append(";");
			
			// rotate back
			s.append("set echo top center;");
			s.append("echo Front view;");
	        s.append("color echo white;");
			s.append("move 0 0 ");
			s.append(360.0-angle);
			s.append(" 0 0 0 0 0 3;");
			
			// show Front view
			s.append("move 0 -90 0 0 0 0 0 0 3;");
			s.append("delay ");
			s.append(delay);
			s.append(";");
			s.append("set echo top center;");
			s.append("echo;");
		} 

		return s.toString();
	}
	
	public String getJmolAnimationOctahedral(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());

	    // show front view
		s.append("set echo top center;");
		s.append("echo Front view: C4-axis;");
        s.append("color echo white;");
        s.append("font echo 24 sanserif;");
        s.append("set echo bottom center;");
		s.append("echo Point group ");
		s.append("O (4,3,2)");
		s.append(";");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo 4-fold rotation;");
        s.append("color echo red;");
		s.append(getJmolSymmetryAxes());
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		
		// rotate around principal axis
		s.append("move 0 0 90 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// align minor axis with z-axis

		// show side view (45 deg.)
		s.append("set echo top center;");
		s.append("echo Diagonal view: C3-axis;");
	    s.append("color echo white;");
		s.append("move 0 45 0 0 0 0 0 0 3;");
		s.append("move 35.3 0 0 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo 3-fold rotation;");
		s.append("color echo royalblue;");
		
//		if (true) {
//			return s.toString();
//		}
		s.append("move 0 0 120 0 0 0 0 0 3;");
		
		// side view
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo Side view: C2-axis;");
	    s.append("color echo white;");
		s.append("move -35.3 0 0 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		
		// 2-fold rotation
		s.append("set echo top center;");
		s.append("echo 2-fold rotation;");
		s.append("color echo royalblue;");
		s.append("move 0 0 180 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 3;");
		s.append("move 0 0 45 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	public String getJmolAnimationIcosahedral(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());

	    // show front view
		s.append("set echo top center;");
		s.append("echo Front view: C5-axis;");
        s.append("color echo white;");
        s.append("font echo 24 sanserif;");
        s.append("set echo bottom center;");
		s.append("echo Point group ");
		s.append("I (5,3,2)");
		s.append(";");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo 5-fold rotation;");
        s.append("color echo red;");
		s.append(getJmolSymmetryAxes());
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		
		// rotate around principal axis
		s.append("move 0 0 72 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// align minor axis with z-axis

		// show side view (45 deg.)
		s.append("set echo top center;");
		s.append("echo Diagonal view: C3-axis;");
	    s.append("color echo white;");
		s.append("move 0 37.5 0 0 0 0 0  0 3;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		
		s.append("set echo top center;");
		s.append("echo 3-fold rotation;");
		s.append("color echo royalblue;");
		s.append("move 0 0 120 0 0 0 0 0 3;");
	
	//	if (true) {
	//		return s.toString();
	//	}
		// side view
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo Side view: C2-axis;");
	    s.append("color echo white;");
		s.append("move 0 0 -120 0 0 0 0 0 1;");
		s.append("move 0 -37.5 0 0 0 0 0 0 1;");
		s.append("move 0 0 72 0 0 0 0 0 1;");
		s.append("move 90 0 0 0 0 0 0 0 1;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		
		// 2-fold rotation
		s.append("set echo top center;");
		s.append("echo 2-fold rotation;");
		s.append("color echo royalblue;");
	    s.append("move 0 0 180 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		
	//	if (true) {
	//		return s.toString();
	//	}

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("move -90 0 0 0 0 0 0 0 3;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	private String getLineJmol(int index, Point3d center, Vector3d axis, String color) {
		StringBuilder s = new StringBuilder();
		s.append("draw l");
		s.append(index);
		s.append(" ");
		s.append("line");
		s.append(" ");

		Vector3d[] vertexes = getPolygonVertices(axis, referenceAxis, center, 2);
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

	private String getPolygonJmol(int index, Point3d center, Vector3d axis, int n, String color) {
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

		Vector3d[] vertexes = getPolygonVertices(axis, referenceAxis, center, n);
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

	private Vector3d[] getPolygonVertices(Vector3d axis, Vector3d referenceAxis, Point3d center, int n) {
		Vector3d perp = new Vector3d(axis);
		// if axis coincides with principal axis, use the reference axis to orient polygon
		if (Math.abs(axis.dot(principalAxis)) > 0.9) {
			perp.set(referenceAxis);
		}
		perp.scale(2);		
		perp.cross(perp, principalAxis);

		AxisAngle4d axisAngle = new AxisAngle4d(axis, 0);
		Vector3d[] vectors = new Vector3d[n];
		Matrix4d m = new Matrix4d();

		for (int i = 0; i < n; i++) {
			//		axisAngle.setAngle(i * 2 * Math.PI/n);
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

	private static double getTrace(Matrix4d matrix) {
		GMatrix m = new GMatrix(4,4);
		m.set(matrix);
		System.out.println("Trace: " + m.trace());
		if (m.trace() <= 0) {
			System.out.println(matrix);
		}
		return m.trace();
	}
}
