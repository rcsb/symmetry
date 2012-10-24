package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.GMatrix;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

public class AxisTransformation {
	private static Vector3d X_AXIS = new Vector3d(1,0,0);
	private static Vector3d NX_AXIS = new Vector3d(-1,0,0);
	private static Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static Vector3d NY_AXIS = new Vector3d(0,-1,0);
	private static Vector3d Z_AXIS = new Vector3d(0,0,1);
	private static Vector3d NZ_AXIS = new Vector3d(0,0,-1);
	
	private static String PRINCIPAL_AXIS_COLOR = "orange";
	private static String MINOR_AXIS_COLOR = "cyan"; // royalblue

	private Subunits subunits = null;
	private RotationGroup rotationGroup = null;

	private Matrix4d transformationMatrix = new Matrix4d();
	private Matrix4d reverseTransformationMatrix = new Matrix4d();
	private Vector3d principalAxis = new Vector3d();
	private Vector3d referenceAxis = new Vector3d();

	private double xMin = Double.MAX_VALUE;
	private double xMax = Double.MIN_VALUE;
	private double yMin = Double.MAX_VALUE;
	private double yMax = Double.MIN_VALUE;
	private double zMin = Double.MAX_VALUE;
	private double zMax = Double.MIN_VALUE;

	boolean modified = true;

	public AxisTransformation(Subunits subunits, RotationGroup rotationGroup) {
		this.subunits = subunits;
		this.rotationGroup = rotationGroup;
	}

	public Matrix4d getTransformation() {
		if (modified) {
			run();   
		}
		return transformationMatrix;
	}
	
	public Matrix4d getReverseTransformation() {
		if (modified) {
			run();   
		}
		return reverseTransformationMatrix;
	}

	private void run () {
		calcTransformation();
		calcBoundaries();
		modified = false;
	}

	private Matrix4d calcTransformation() {
		if (subunits == null || subunits.getSubunitCount() == 0) {
			// if the structure doesn't contain protein subunits, i.e.,
			// a DNA structure, return the identity matrix.
			transformationMatrix.setIdentity();
		} else if (rotationGroup.getPointGroup().equals("C1")) {
			transformationMatrix = getTransformationByInertiaAxes();
		} else {
		//	transformationMatrix = getTransformationByInertiaAxes();
			transformationMatrix = getTransformationBySymmetryAxes();
			// for cyclic geometry, set a canonical view for the Z direction
			if (rotationGroup.getPointGroup().startsWith("C")) {
				setZDirection(transformationMatrix);
			}
		}

		return transformationMatrix;
	}

	private Matrix4d getTransformationBySymmetryAxes() {
		Point3d[] refPoints = new Point3d[6];

		// second reference point is along the principal rotation axis
		principalAxis = getPrincipalRotationAxis();
		principalAxis.normalize();

		// superposition cannot handle 180 degree
		if (Z_AXIS.dot(principalAxis) < -0.9) {
			principalAxis.negate();
		}
		// third reference point is along an orthogonal rotation axis (if available), otherwise,
		// an orthogonal vector is used.
		// if rotation group has orthogonal axis, use it for alignment
		referenceAxis = getMinorRotationAxis();
		if (referenceAxis == null) {
			referenceAxis = getSubunitReferenceAxisNew();
		}
		System.out.println("z.dot.p: " + Z_AXIS.dot(principalAxis));
		System.out.println("ref.dot.principal" + Math.abs(referenceAxis.dot(principalAxis)));
		System.out.println("ref-principal angle" + Math.toDegrees(referenceAxis.angle(principalAxis)));
		if (Math.abs(referenceAxis.dot(principalAxis)) > 0.000001) {
			// TODO new reference axis should point into the same direction as original reference axis
			referenceAxis.cross(principalAxis, referenceAxis); // make it perpendicular
			referenceAxis.cross(referenceAxis, principalAxis);
			System.out.println("non-perpendicular axis");
			System.out.println("ref-principal angle corrected:" + Math.toDegrees(referenceAxis.angle(principalAxis)));
		}
		referenceAxis.normalize();
		if (Y_AXIS.dot(referenceAxis) < -0.9) {
			referenceAxis.negate();
		}
		Vector3d perpendicularAxis = new Vector3d();
		perpendicularAxis.cross(principalAxis, referenceAxis);
		
		refPoints[0] = new Point3d(principalAxis);
		refPoints[1] = new Point3d(refPoints[0]);
		refPoints[1].negate();

		refPoints[2] = new Point3d(referenceAxis);
		refPoints[3] = new Point3d(refPoints[2]);
		refPoints[3].negate();
		
		refPoints[4] = new Point3d(perpendicularAxis);
		refPoints[5] = new Point3d(refPoints[4]);
		refPoints[5].negate();

//		System.out.println("Reference points");
//		for (Point3d p:refPoints) {
//			System.out.println(p);
//		}
		//  y,z axis centered at the centroid of the subunits
		Point3d[] coordPoints = new Point3d[6];
		coordPoints[0] = new Point3d(Z_AXIS);
		coordPoints[1] = new Point3d(NZ_AXIS);
		coordPoints[2] = new Point3d(Y_AXIS);
		coordPoints[3] = new Point3d(NY_AXIS);
		coordPoints[4] = new Point3d(NX_AXIS);
		coordPoints[5] = new Point3d(X_AXIS);
//		System.out.println("coord points");
//		for (Point3d p:coordPoints) {
//			System.out.println(p);
//		}
 
		// align principal axis with z axis and perpendicular axis with the y axis
        Point3d[] ref1 = new Point3d[refPoints.length];
        for (int i = 0; i < ref1.length; i++) {
        	ref1[i] = new Point3d(refPoints[i]);
        }
        long t1 = System.nanoTime();
		Matrix4d matrix = SuperPosition.superposeAtOrigin(refPoints, coordPoints);
        long t2 = System.nanoTime();
        System.out.println("old superposition: " + (t2-t1));
        System.out.println("old matrix: " + matrix);
        long t3 = System.nanoTime();
		SuperPosition s = new SuperPosition(false, 0);
	    s.calcQCPSuperpositionRotationOnly(ref1, coordPoints);
		Matrix4d m3 = s.getTransformationMatrix();
		long t4 = System.nanoTime();
        System.out.println("new superposition: " + (t4-t3));
        System.out.println("new matrix: " + m3);
       
		calcReverseMatrix(matrix);

		return matrix;
	}

	private void calcReverseMatrix(Matrix4d matrix) {
		reverseTransformationMatrix.set(matrix);
		// reset translation vector
		reverseTransformationMatrix.setColumn(3, 0, 0, 0, 1);
		reverseTransformationMatrix.transpose();
		// set translation vector to centroid
		Point3d centroid = subunits.getCentroid();
		reverseTransformationMatrix.setColumn(3, centroid.x, centroid.y, centroid.z, 1);
		System.out.println("ReverseTransformation: " + reverseTransformationMatrix);
		System.out.println("Centroid: " + subunits.getCentroid());
	}

	/**
	 * Returns vector from largest subunit to centroid of complex
	 * TODO would it be better to use an inertia axis instead if it is orthogonal to the
	 * principal rotation axis??
	 * @return
	 */	
	private Vector3d getSubunitReferenceAxisNew() {
		Vector3d[] inertiaVectors = calcPrincipalAxes();
		for (Vector3d v: inertiaVectors) {
			if (Math.abs(principalAxis.dot(v)) < 0.1) {
				return v;
			}
		}
		return inertiaVectors[1];
	}

	/**
	 * Returns the z-radius
	 * @return double radius in z-direction
	 */
	private double getZRadius() {
		System.out.println("zMax: " + zMax + " zMin" + zMin);
		return 0.5 * (zMax - zMin); // half of dimension along z-axis (principal rotation axis)
	}
	
	/**
	 * Returns the radius for drawing the minor rotation axis in the xy-plane
	 * @return double radius in xy-plane
	 */
	private double getXYRadius() {
//	    double r1 = 0.25 * Math.max((xMax-xMin), (yMax-yMin)); // max radius for rotation in xy plane
//	    double r2 = 0.25 * Math.sqrt((xMax-xMin)*(xMax-xMin) + (yMax-yMin)*(yMax-yMin));
	//	return r1 + r2;
		return  0.5 * Math.max((xMax-xMin), (yMax-yMin));
	}
	
	/**
	 * Returns the length of side of an n-fold regular polygon given its inradius
	 * @return
	 */
	private static double polygonSide(int n, double inradius) {
	    return inradius * 2 * Math.tan(Math.PI/n);
	}
	
	/**
	 * Returns the radius of an n-fold regular polygon given the length of a side
	 * @return
	 */
	private static double polygonRadius(int n, double length) {
	    return length / (2 * Math.sin(Math.PI/n));
	}

	/**
	 * Returns a transformation matrix for prisms to draw for Cn structures.
	 * The center in this matrix is the geometric center, rather then the centroid
	 * In Cn structures those are usually not the same.
	 * @return
	 */
	private Matrix4d calcPrismTransformation() {
		Matrix4d geometricCentered = new Matrix4d(reverseTransformationMatrix);
		// reset translation vector
		geometricCentered.setColumn(3, 0, 0, 0, 1);
		// set translation vector to centroid
		Point3d geometricCenter = calcGeometricCenter();
		
		geometricCentered.setColumn(3, geometricCenter.x, geometricCenter.y, geometricCenter.z, 1);
		System.out.println("GeometricCentered: " + geometricCentered);
		System.out.println("Centroid: " + subunits.getCentroid());
		return geometricCentered;
	}

	private Point3d calcGeometricCenter() {
		Point3d geometricCenter = new Point3d();
	
		if (rotationGroup.getPointGroup().startsWith("C")) {
	//		geometricCenter = new Point3d(xMin + (xMax-xMin)*0.5, yMin +(yMax-yMin)*0.5, zMin + (zMax-zMin)*0.5);
	//		geometricCenter = new Point3d(0, 0, zMin + (zMax-zMin)*0.5);
	//		transformationMatrix.transform(geometricCenter);
		}
		
		geometricCenter.add(subunits.getCentroid());
		return geometricCenter;
	}
 	/**
	 * Calculates the min and max boundaries of the structure after it has been
	 * transformed into its canonical orientation.
	 */
	private void calcBoundaries() {
		if (subunits == null || subunits.getSubunitCount() == 0) {
			return;
		}
		xMin = Double.MAX_VALUE;
		xMax = Double.MIN_VALUE;
		yMin = Double.MAX_VALUE;
		yMax = Double.MIN_VALUE;
		zMin = Double.MAX_VALUE;
		zMax = Double.MIN_VALUE;

		Point3d probe = new Point3d();
		Point3d centroid = subunits.getCentroid();

		for (Point3d[] list: subunits.getTraces()) {
			for (Point3d p: list) {
				probe.sub(p, centroid); // apply transformation at origin of coordinate system
				transformationMatrix.transform(probe);

				xMin = Math.min(xMin, probe.x);
				xMax = Math.max(xMax, probe.x);
				yMin = Math.min(yMin, probe.y);
				yMax = Math.max(yMax, probe.y);
				zMin = Math.min(zMin, probe.z);
				zMax = Math.max(zMax, probe.z);
			}
		}
	}

	/*
	 * Modifies the rotation part of the transformation axis for
	 * a Cn symmetric complex, so that the narrower end faces the
	 * viewer, and the wider end faces away from the viewer. Example: 3LSV
	 */
	private Matrix4d setZDirection(Matrix4d matrix) {
		calcBoundaries();

		// TODO can we just use the z-min/z-max distances to create canonical view?
		Point3d probe = new Point3d();
		Point3d centroid = subunits.getCentroid();
		// TODO centroid required?
		centroid = new Point3d();

		// find the center along the z-axis
		double center = zMin + (zMax - zMin)/2;
		double sum1 = 0;
		double sum2 = 0;
		for (Point3d[] list: subunits.getTraces()) { // loop over all subunits
			for (Point3d p: list) {			
				// align points with z-axis (principal rotation axis)
				probe.sub(p, centroid);
				transformationMatrix.transform(probe);
				
				// calculate the distance square for each
				// point from the z-axis. Sum up the distance 
				// squares for front and back half of the structure.
				if (probe.z < center) {
					sum1 += probe.x*probe.x + probe.y*probe.y;
				} else {
					sum2 += probe.x*probe.x + probe.y*probe.y;
				}
			}
		}

		// if the wider part (large sum) faces the front, apply
		// a 180 degree rotation around the y-axis so that the 
		// narrower part faces the viewer. This new orientation 
		// provides the least occluded view along the z-axis.
		if (sum2 > sum1) {
			Matrix4d rot = new Matrix4d();
			rot.set(new AxisAngle4d(0, 1, 0, -Math.PI));
			matrix.mul(rot);
			System.out.println("Changing z direction");
		}

		return matrix;
	}

	/**
	 * Returns a normalized vector that represents the
	 * principal rotation axis.
	 * @return principal rotation axis
	 */
	private Vector3d getPrincipalRotationAxis() {
		Rotation rotation = rotationGroup.getRotation(0); // the rotation around the principal axis is the first rotation
		AxisAngle4d axisAngle = rotation.getAxisAngle();
		Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		v.normalize();
		return v;
	}

	/**
	 * Returns a normalized vector that represents the
	 * minor rotation axis, or null if there is no minor rotation axis in
	 * the case of cyclic symmetry.
	 * @return minor rotation axis, return null if minor rotation axis does not exists
	 */
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

	/**
	 * Returns the fold of the highest order rotation axis
	 * @return highest order rotation
	 */
	private int getPrincipalRotationAxisFold() {
		if (rotationGroup.getOrder() == 0) {
			return 1;
		}
		Rotation rotation = rotationGroup.getRotation(0); // the rotation around the principal axis is the first rotation
		return rotation.getFold();
	}

	/**
	 * Returns the fold of a minor order rotation axis
	 * @return minor order rotation
	 */
	private int getMinorRotationAxisFold() {
		// find axis that is not the rotation principal axis (direction = 1)
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
			if (rotationGroup.getRotation(i).getDirection() == 1) {
				return rotationGroup.getRotation(i).getFold();
			}
		}
		return 1;
	}

	private Vector3d[] calcPrincipalAxes() {
		MomentsOfInertia moi = new MomentsOfInertia();

		for (Point3d[] list: subunits.getTraces()) {
			for (Point3d p: list) {
				moi.addPoint(p, 1.0);
			}
		}
		return moi.getPrincipalAxes();
	}

	private Matrix4d getTransformationByInertiaAxes() {
		Vector3d[] inertiaVectors = calcPrincipalAxes();
		drawInertiaAxes(subunits.getCentroid(), inertiaVectors);

		if (Y_AXIS.dot(inertiaVectors[0]) < -0.9) {
			inertiaVectors[0].negate();
		}
		if (X_AXIS.dot(inertiaVectors[1]) < -0.9) {
			inertiaVectors[1].negate();
		}
		// orientation of subunit inertia vectors center at the centroid of the subunits
		Point3d[] refPoints = new Point3d[inertiaVectors.length*2];
		refPoints[0] = new Point3d(inertiaVectors[0]);
		refPoints[1] = new Point3d(refPoints[0]);
		refPoints[1].negate();
		
		refPoints[2] = new Point3d(inertiaVectors[1]);
		refPoints[3] = new Point3d(refPoints[2]);
		refPoints[3].negate();
		
		refPoints[4] = new Point3d(inertiaVectors[2]);
		refPoints[5] = new Point3d(refPoints[4]);
		refPoints[5].negate();
		

		// x,y,z axis center at the centroid of the subunits
		Point3d[] coordPoints = new Point3d[6];
		coordPoints[0] = new Point3d(Y_AXIS);
		coordPoints[1] = new Point3d(NY_AXIS);
		coordPoints[2] = new Point3d(X_AXIS);
		coordPoints[3] = new Point3d(NX_AXIS);
		coordPoints[4] = new Point3d(NZ_AXIS);
		coordPoints[5] = new Point3d(Z_AXIS);
		
		// align inertia axis with x,y,z axis
		Matrix4d matrix = SuperPosition.superposeAtOrigin(refPoints, coordPoints);
		
		calcReverseMatrix(matrix);
		return matrix;
	}

	private void drawInertiaAxes(Point3d center, Vector3d[] axes) {
		StringBuilder s = new StringBuilder();

		for (int i = 0; i < axes.length; i++) {
			s.append("draw l");
			s.append(i);
			s.append(" ");
			s.append("line");
			s.append(" ");
			s.append("{");
			s.append(jMolFloat(center.x));
			s.append(" ");
			s.append(jMolFloat(center.y));
			s.append(" ");
			s.append(jMolFloat(center.z));
			s.append("}");
			s.append("{");
			Vector3d v = new Vector3d(axes[i]);
			v.scale(20);
			v.add(center);
			s.append(jMolFloat(v.x));
			s.append(" ");
			s.append(jMolFloat(v.y));
			s.append(" ");
			s.append(jMolFloat(v.z));
			s.append("}");
			s.append(" width 0.5 ");
			s.append(" color white");
			s.append(";");
		}
        System.out.println(s);
	}

	/**
	 * Returns a Jmol script to rotate a structure into
	 * a canonical orientation. For symmetric structures,
	 * the highest order axis is aligned with the z-axis, 
	 * and for asymmetric structures, the principal axis of
	 * inertia is aligned with the y-axis.
	 * @return Jmol script
	 */
	public String getJmolTransformation() {
		if (modified) {
			run();
		}
		Quat4d q = new Quat4d();
		transformationMatrix.get(q);
		Point3d centroid = subunits.getCentroid();
		return "rotate quaternion {" + jMolFloat(q.x) + " " + jMolFloat(q.y) + " " + jMolFloat(q.z) + " " + jMolFloat(q.w) + "}; zoom 80; center {"+ centroid.x + " " + centroid.y + " " + centroid.z + "};";
	}

	public Matrix4d getInverseTransformation() {
		if (modified) {
			run();
		}
		Matrix4d m = new Matrix4d(transformationMatrix);
		m.invert();
		return m;
	}
	
	public String getJmolSymmetryAxes() {
		StringBuilder s = new StringBuilder();

		int n = rotationGroup.getOrder();
		if (n == 0) {
			return s.toString();
		}
        int maxFold = getPrincipalRotationAxisFold();

		float diameter = 0.5f;
		double radius = 0;
		String color = "";

		for (int i = 0; i < n; i++) {
			Rotation rotation = rotationGroup.getRotation(i);
			int direction =  rotation.getDirection();

			// don't draw redundant n-fold rotations around principal axis
			if (i > 0 && direction == 0) {
				continue;
			}

			if (direction == 0) {
				radius = getZRadius()*1.5; // principal axis uses z-dimension
				color = PRINCIPAL_AXIS_COLOR;
				diameter = 0.25f;
			} else {
				// convert inradius into length of side of polygon
				if (maxFold > 2) {
					double side = polygonSide(maxFold, getXYRadius());
					// calculate radius of polygon from side of polygon
					radius = 1.2 * polygonRadius(maxFold, side);
					
				} else {
					radius = 1.2 * getXYRadius();
				}
				color = MINOR_AXIS_COLOR;
				diameter = 0.1f;
			}
			
			Point3d center = calcGeometricCenter();
			AxisAngle4d axisAngle = rotation.getAxisAngle();
			Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
			
			s.append(JmolGeometricObjects.getSymmetryAxis(i, rotationGroup.getPointGroup(), rotation.getFold(), principalAxis, referenceAxis, radius, diameter, color, center, axis));
						
			if (color.isEmpty()) {
			Point3d p1 = new Point3d(axis);
			p1.scaleAdd(-radius, calcGeometricCenter());

			Point3d p2 = new Point3d(axis);
			p2.scaleAdd(radius, calcGeometricCenter());

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
			p1.scaleAdd(-radius*0.95, calcGeometricCenter());

			p2 = new Point3d(axis);
			p2.scaleAdd(radius*0.95, calcGeometricCenter());

			if (rotation.getFold() == 2) {
				s.append(getLineJmol(i, p1, axis, color));
				s.append(getLineJmol(i + n, p2, axis, color));
			} else {
				double polygonRadius = 2;
				s.append(getPolygonJmol(i, p1, axis, rotation.getFold(), color, polygonRadius));
				s.append(getPolygonJmol(i + n, p2, axis, rotation.getFold(), color, polygonRadius));
			}
			}
		}
		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure in Jmol
	 * into canonical orientations and animates the rotation
	 * around unique symmetry axes.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	public String getJmolAnimation(int delay) {
		String animation = "";
		String pointGroup = rotationGroup.getPointGroup();
		if (pointGroup.startsWith("C")) {
			animation = getJmolAnimationCyclic(delay);
		} else if (pointGroup.startsWith("D")) {
			animation = getJmolAnimationDihedral(delay);
		} else if (pointGroup.equals("T")) {
			animation = getJmolAnimationTetrahedral(delay);
		} else if (pointGroup.equals("O")) {
			animation = getJmolAnimationOctahedral(delay);
		} else if (pointGroup.equals("I")) {
			animation = getJmolAnimationIcosahedral(delay);
		}
		return animation;
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * cyclic (Cn) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the Cn axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
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
		// TODO could use the inertia axes to do rotations for C1 symmetry
		if (fold > 1) {
			// show symmetry axes
			s.append("set echo top center;");
			s.append("echo ");
			s.append(fold);
			s.append("-fold rotation;");
			s.append("color echo ");
			s.append(PRINCIPAL_AXIS_COLOR);
			s.append(";");
			s.append(getJmolSymmetryAxes());
			double radius = 0;
			if (fold > 2) {
				double side = polygonSide(fold, getXYRadius());
				// calculate radius of polygon from side of polygon
				radius = polygonRadius(fold, side);
				System.out.println("r(xy): " + getXYRadius() + "r(p): " + radius);
			} else {
				radius = getXYRadius();
			}
			// TODO the center for Cn should be the geometric center, not the centroid (part of revereseTransformationMatrix)
			s.append(JmolGeometricObjects.getJmolWirePrism(100, fold, getZRadius()*2, radius, calcPrismTransformation()));
			s.append("delay ");
			s.append(delay); 
			s.append(";");

			// rotate around principal axis
			fold  = getPrincipalRotationAxisFold();
			float angle = 360.0f/fold;
			s.append("move 0 0 ");
			s.append(angle);
			s.append(" 0 0 0 0 0 4;");
			s.append("delay ");
			s.append(delay);
			s.append(";");
		}

		// show Side view 1
		s.append("set echo top center;");
		s.append("echo Side view 1;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Back view
		s.append("set echo top center;");
		s.append("echo Back view;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Side view 2
		s.append("set echo top center;");
		s.append("echo Side view 2;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * dihedral (Dn) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the Cn and perpendicular
	 * C2 axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
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
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
		s.append(getJmolSymmetryAxes());
		double radius = 0;
		if (fold > 2) {
			double side = polygonSide(fold, getXYRadius());
			// calculate radius of polygon from side of polygon
			radius = polygonRadius(fold, side);
			System.out.println("r(xy): " + getXYRadius() + "r(p): " + radius);
		} else {
			radius = getXYRadius();
		}
		s.append(JmolGeometricObjects.getJmolWirePrism(100, fold, getZRadius()*2, radius, reverseTransformationMatrix));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		float angle = 360.0f/fold;
		s.append("move 0 0 ");
		s.append(angle);
		s.append(" 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show side view
		s.append("set echo top center;");
		s.append("echo Side view;");
		s.append("color echo white;");
		// align minor axis with z-axis
		s.append("move 90 0 0 0 0 0 0 0 4;");

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
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 ");
		s.append(angle);
		s.append(" 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// rotate back
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("move 0 0 ");
		s.append(360.0-angle);
		s.append(" 0 0 0 0 0 0;");

		// show Front view
		s.append("move -90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * tetrahedral (T) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the C3 and C2 axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	public String getJmolAnimationTetrahedral(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());

		// show front view
		s.append("set echo top center;");
		s.append("echo Front view: C3-axis;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("set echo bottom center;");
		s.append("echo Point group ");
		s.append("T (2,3)");
		s.append(";");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo 3-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
		s.append(getJmolSymmetryAxes());

		double radius = Math.max(getZRadius(), getXYRadius());
		System.out.println("zMaxRad: " + getZRadius());
		System.out.println("xyMaxRad: " + getXYRadius());
		radius = JmolGeometricObjects.tetrahedronInRadiusToOutRadius(radius);
		s.append(JmolGeometricObjects.getJmolWireTetrahedron(100, radius, reverseTransformationMatrix));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		s.append("move 0 0 120 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		s.append("set echo top center;");
		s.append("echo Edge view: C2-axis;");
		s.append("color echo white;");
		double tetrahedralAngle = Math.toDegrees(Math.acos(-1.0/3.0));
		double angle = 180 - 0.5 * tetrahedralAngle;
		System.out.println("tetaangle: " + angle);
		s.append("move ");
		s.append(jMolFloat(angle));
		s.append(" 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// rotate around 2-fold axis
		s.append("set echo top center;");
		s.append("echo 2-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 180 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// rotate back to front view
		s.append("echo Front view: C3-axis;");
		s.append("color echo white;");
		s.append("move 0 0 -180 0 0 0 0 0 0;");
		s.append("move ");
		s.append(jMolFloat(-angle));
		s.append(" 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * octahedral (O) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the C4, C3, and C2 axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
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
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
		s.append(getJmolSymmetryAxes());
		double radius = Math.max(getZRadius(), getXYRadius());
		radius = JmolGeometricObjects.octahedronInRadiusToOutRadius(radius);
		s.append(JmolGeometricObjects.getJmolWireOctahedron(100, radius, reverseTransformationMatrix));
		
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		s.append("move 0 0 90 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show side view (45 deg.)
		s.append("set echo top center;");
		s.append("echo Diagonal view: C3-axis;");
		s.append("color echo white;");
		s.append("move 0 45 0 0 0 0 0 0 0;");
		s.append("move 35.3 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo 3-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 120 0 0 0 0 0 4;");

		// side view
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo Side view: C2-axis;");
		s.append("color echo white;");
		s.append("move -35.3 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// 2-fold rotation
		s.append("set echo top center;");
		s.append("echo 2-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 180 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 0;");
		s.append("move 0 0 45 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * icosahedral (I) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the C5, C3, and C2 axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
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
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
		s.append(getJmolSymmetryAxes());
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		s.append("move 0 0 72 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// align minor axis with z-axis

		// show side view (45 deg.)
		s.append("set echo top center;");
		s.append("echo Diagonal view: C3-axis;");
		s.append("color echo white;");
		s.append("move 0 37.5 0 0 0 0 0  0 4;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		s.append("set echo top center;");
		s.append("echo 3-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 120 0 0 0 0 0 4;");

		// side view
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo Side view: C2-axis;");
		s.append("color echo white;");
		s.append("move 0 0 -120 0 0 0 0 0 0;");
		s.append("move 0 -37.5 0 0 0 0 0 0 0;");
		s.append("move 0 0 72 0 0 0 0 0 0;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// 2-fold rotation
		s.append("set echo top center;");
		s.append("echo 2-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 180 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("move -90 0 0 0 0 0 0 0 4;");
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

		Vector3d[] vertexes = getPolygonVertices(axis, referenceAxis, center, 2, 2);
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

	private String getPolygonJmol(int index, Point3d center, Vector3d axis, int n, String color, double radius) {
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

		Vector3d[] vertexes = getPolygonVertices(axis, referenceAxis, center, n, radius);
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

	private Vector3d[] getPolygonVertices(Vector3d axis, Vector3d referenceAxis, Point3d center, int n, double radius) {
		Vector3d perp = new Vector3d(axis);
//		System.out.println("AxisTransformation: radius: " + radius);
//		System.out.println("AxisTransformation: center: " + center);
//		System.out.println("AxisTransformation: principal axis: " + principalAxis);
//		System.out.println("AxisTransformation: reference axis: " + referenceAxis);
//		System.out.println("AxisTransformation: axis: " + axis);
		// if axis coincides with principal axis, use the reference axis to orient polygon
		// TODO need to check alignment with y-axis for all cases of symmetry
		if (Math.abs(axis.dot(principalAxis)) > 0.9) {
			perp.set(referenceAxis);
		} else {
			perp.cross(perp, principalAxis);
		}
//		System.out.println("AxisTransformation: perp. axis: " + referenceAxis);
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
