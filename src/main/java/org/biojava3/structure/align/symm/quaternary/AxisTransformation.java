package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

public class AxisTransformation {
	private static Vector3d X_AXIS = new Vector3d(1,0,0);
	private static Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static Vector3d Z_AXIS = new Vector3d(0,0,1);
	
	private static String PRINCIPAL_AXIS_COLOR = "darkorange";
	private static String MINOR_AXIS_COLOR = "darkturquoise"; // royalblue

	private Subunits subunits = null;
	private RotationGroup rotationGroup = null;

	private Matrix4d transformationMatrix = new Matrix4d();
	private Matrix4d reverseTransformationMatrix = new Matrix4d();
	private Vector3d principalAxis = new Vector3d();
	private Vector3d referenceAxis = new Vector3d();
	private Vector3d[] principalAxes = null;

	private double xMin = Double.MAX_VALUE;
	private double xMax = Double.MIN_VALUE;
	private double yMin = Double.MAX_VALUE;
	private double yMax = Double.MIN_VALUE;
	private double zMin = Double.MAX_VALUE;
	private double zMax = Double.MIN_VALUE;
	private double xyRadiusMax = Double.MIN_VALUE;

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
		} else if (rotationGroup.getPointGroup().equals("C2")) {
			transformationMatrix = getTransformationByC2SymmetryAxes();;
		} else {
			transformationMatrix = getTransformationBySymmetryAxes();
			// for cyclic geometry, set a canonical view for the Z direction
			if (rotationGroup.getPointGroup().startsWith("C")) {
				setZDirection(transformationMatrix);
			}
		}

		return transformationMatrix;
	}

	private Matrix4d getTransformationBySymmetryAxes() {
		principalAxis = getPrincipalRotationAxis();
		referenceAxis = getMinorRotationAxis();
		if (referenceAxis == null) {
			referenceAxis = getSubunitReferenceAxisNew();
		}
		// make sure reference axis is perpendicular
		referenceAxis.cross(principalAxis, referenceAxis);
		referenceAxis.cross(referenceAxis, principalAxis);

		Point3d[] refPoints = new Point3d[2];
		refPoints[0] = new Point3d(principalAxis);
		refPoints[1] = new Point3d(referenceAxis);


		//  y,z axis centered at the centroid of the subunits
		Point3d[] coordPoints = new Point3d[2];
		coordPoints[0] = new Point3d(Z_AXIS);
		coordPoints[1] = new Point3d(Y_AXIS);

		Matrix4d matrix = alignAxes(refPoints, coordPoints);

		calcReverseMatrix(matrix);

		return matrix;
	}
	
	private Matrix4d getTransformationByC2SymmetryAxes() {
		principalAxis = getPrincipalRotationAxis();
		referenceAxis = getMinorRotationAxis();
		if (referenceAxis == null) {
			referenceAxis = getSubunitReferenceAxisNew();
		}
		
		// make sure reference axis is perpendicular
		referenceAxis.cross(principalAxis, referenceAxis);
		referenceAxis.cross(referenceAxis, principalAxis);
		
		Point3d[] refPoints = new Point3d[2];
		refPoints[0] = new Point3d(principalAxis);
		refPoints[1] = new Point3d(referenceAxis);
		
		Point3d[] coordPoints = new Point3d[2];
		coordPoints[0] = new Point3d(Y_AXIS);
		coordPoints[1] = new Point3d(X_AXIS);
 
    	Matrix4d matrix = alignAxes(refPoints, coordPoints);  
		calcReverseMatrix(matrix);

		return matrix;
	}
	
	private Matrix4d alignAxes(Point3d[] refPoints, Point3d[] coordPoints) {
//		System.out.print("P0: " + refPoints[0] + " ");
	
		Vector3d v1 = new Vector3d(refPoints[0]);
//		System.out.println("V0 len: " + v1.length());
		Vector3d v2 = new Vector3d(coordPoints[0]);
		Matrix4d m1 = new Matrix4d();
		AxisAngle4d a = new AxisAngle4d();
		Vector3d axis = new Vector3d();
		axis.cross(v1,v2);
		a.set(axis, v1.angle(v2));
		m1.set(a);
		for (int i = 0; i < 2; i++) {
			Point3d v = refPoints[i];
			m1.transform(v);
		}
		System.out.println(refPoints[0]);
		v1 = new Vector3d(refPoints[1]);
//		System.out.println("V1 len: " + v1.length());
		v2 = new Vector3d(coordPoints[1]);
//		System.out.print("P1: " + refPoints[1] + " ");
		axis.cross(v1,v2);
		a.set(axis, v1.angle(v2));
		Matrix4d m2 = new Matrix4d();
		m2.set(a);
		for (Point3d v: refPoints) {
			m2.transform(v);
		}
		System.out.println(refPoints[1]);
		if (SuperPosition.rmsd(refPoints, coordPoints) > 0.01) {
			System.out.println("Warning: aligment with coordiante system is off. RMSD: " + SuperPosition.rmsd(refPoints, coordPoints));
		}
		m2.mul(m1);
		return m2;
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
//		System.out.println("Centroid: " + subunits.getCentroid());
	}

	/**
	 * Returns vector from largest subunit to centroid of complex
	 * TODO would it be better to use an inertia axis instead if it is orthogonal to the
	 * principal rotation axis??
	 * @return
	 */	
	private Vector3d getSubunitReferenceAxisNew() {
		if (rotationGroup.getPointGroup().equals("C2")) {
			Vector3d vr = new Vector3d(subunits.getCentroid());
			vr.sub(subunits.getOriginalCenters().get(0));
			vr.normalize();
			return vr;
		}		
		
		Vector3d[] inertiaVectors = calcPrincipalAxes();
		Vector3d vmin = null;
		double dotMin = 1.0;
		for (Vector3d v: inertiaVectors) {
			if (Math.abs(principalAxis.dot(v)) < dotMin) {
				dotMin = Math.abs(principalAxis.dot(v));
				vmin = new Vector3d(v);
			}
		}
		return vmin;
	}

	/**
	 * Returns the x-radius
	 * @return double radius in z-direction
	 */
	private double getXRadius() {
		return 0.5 * (xMax - xMin); // half of dimension along z-axis (principal rotation axis)
	}
	/**
	 * Returns the y-radius
	 * @return double radius in z-direction
	 */
	private double getYRadius() {
		return 0.5 * (yMax - yMin); // half of dimension along z-axis (principal rotation axis)
	}
	
	/**
	 * Returns the z-radius
	 * @return double radius in z-direction
	 */
	private double getZRadius() {
		return 0.5 * (zMax - zMin); // half of dimension along z-axis (principal rotation axis)
	}
	
	/**
	 * Returns the radius for drawing the minor rotation axis in the xy-plane
	 * @return double radius in xy-plane
	 */
	private double getXYRadius() {
		return xyRadiusMax;
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
	
		if (rotationGroup.getPointGroup().equals("C1")) {
			geometricCenter = new Point3d(xMin + getXRadius(), yMin + getYRadius(), zMin + getZRadius());
			// TODO does this transformation include the translational component??
			transformationMatrix.transform(geometricCenter);
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
				xyRadiusMax = Math.max(xyRadiusMax, Math.sqrt(probe.x*probe.x + probe.y * probe.y));
			}
		}
		System.out.println("xMin: " + xMin + " xMax: " + xMax);
		System.out.println("yMin: " + yMin + " yMax: " + yMax);
		System.out.println("zMin: " + zMin + " zMax: " + zMax);
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
    	if (rotationGroup.getPointGroup().equals("I")) {
			return getReferenceAxisIcosahedral();
		}
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
	
	private Vector3d getReferenceAxisIcosahedral() {
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				double d = v.dot(principalAxis);
				if (rotationGroup.getRotation(i).getFold() == 5) {
					// the dot product of 0.447.. is between to adjacent 5-fold axis
					if (d > 0.447 && d < 0.448) {
						return v;
					}
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

		Point3d[] refPoints = new Point3d[2];
		refPoints[0] = new Point3d(inertiaVectors[0]);
		refPoints[1] = new Point3d(inertiaVectors[1]);

		Point3d[] coordPoints = new Point3d[2];
		coordPoints[0] = new Point3d(Y_AXIS);
		coordPoints[1] = new Point3d(X_AXIS);
		
		// align inertia axis with y-x axis
		Matrix4d matrix = alignAxes(refPoints, coordPoints);
		if (SuperPosition.rmsd(refPoints, coordPoints) > 0.01) {
			System.out.println("Warning: aligment with coordiante system is off. RMSD: " + SuperPosition.rmsd(refPoints, coordPoints));
		}
		
		calcReverseMatrix(matrix);
		return matrix;
	}

	private String jmolDrawInertiaAxes() {
		StringBuilder s = new StringBuilder();
		Point3d centroid = calcGeometricCenter();
		if (principalAxes == null) {
			principalAxes = calcPrincipalAxes();
		}

		for (int i = 0; i < principalAxes.length; i++) {
			s.append("draw l");
			s.append(200+i);
			s.append(" ");
			s.append("line");
			s.append(" ");
			Point3d v1 = new Point3d(principalAxes[i]);
			if (i == 0) {
				v1.scale(1.2*getYRadius());
			} else if (i == 1) {
				v1.scale(1.2*getXRadius());
			} else if (i == 2) {
				v1.scale(1.2*getZRadius());
			}
			Point3d v2 = new Point3d(v1);
			v2.negate();
			v1.add(centroid);
			v2.add(centroid);
			s.append(getJmolPoint(v1));
			s.append(getJmolPoint(v2));
			s.append(" width 0.5 ");
			s.append(" color white");
			s.append(";");
		}
        return s.toString();
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
		if (rotationGroup.getPointGroup().equals("C1")) {
			centroid = calcGeometricCenter();
		}
		int zoom = 80;
		if (rotationGroup.getPointGroup().equals("T")) {
			zoom = 60;
		} else if (rotationGroup.getPointGroup().equals("O") || rotationGroup.getPointGroup().equals("I")) {
			zoom = 80;
		}
		return "rotate quaternion {" + jMolFloat(q.x) + " " + jMolFloat(q.y) + " " + jMolFloat(q.z) + " " + jMolFloat(q.w) + "}; zoom " + zoom + "; center {"+ centroid.x + " " + centroid.y + " " + centroid.z + "};";
	//	return "rotate quaternion {" + jMolFloat(q.x) + " " + jMolFloat(q.y) + " " + jMolFloat(q.z) + " " + jMolFloat(q.w) + "}; zoom " + zoom + ";";

	}

//	public Matrix4d getInverseTransformation() {
//		if (modified) {
//			run();
//		}
//		Matrix4d m = new Matrix4d(transformationMatrix);
//		m.invert();
//		return m;
//	}
	
	public String getJmolSymmetryAxes() {
		StringBuilder s = new StringBuilder();

		int n = rotationGroup.getOrder();
		if (n == 0) {
			return s.toString();
		}

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

			if (rotationGroup.getPointGroup().startsWith("C") || rotationGroup.getPointGroup().startsWith("D")) {
				if (direction == 0) {
					radius = 1.2 *getZRadius(); // principal axis uses z-dimension
				} else {
					radius = 1.1 * getXYRadius();
				} 
			} else if (rotationGroup.getPointGroup().equals("T")) {
				radius = 0.9 * Math.max(getZRadius(), getXYRadius());
			} else { 
				radius = Math.max(getZRadius(), getXYRadius()); // for O, I point group
			}

			if (direction == 0) { // for principal axis
	//			radius = 1.2 * radius;
				color = PRINCIPAL_AXIS_COLOR;
				diameter = 0.5f;
			} else { // for all other axes
	//			radius = 1.1 * radius;
				color = MINOR_AXIS_COLOR;
				diameter = 0.25f;
			} 

			Point3d center = calcGeometricCenter();
			AxisAngle4d axisAngle = rotation.getAxisAngle();
			Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);

			s.append(JmolGeometricObjects.getSymmetryAxis(i*100, direction, rotationGroup.getPointGroup(), rotation.getFold(), principalAxis, referenceAxis, radius, diameter, color, center, axis));

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
		if (pointGroup.equals("C1")) {
			animation = getJmolAnimationC1(delay);
		} else if (pointGroup.startsWith("C2")) {
			animation = getJmolAnimationC2(delay);
		//	animation = getJmolAnimationCyclic(delay);
		} else if (pointGroup.startsWith("C")) {
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
	 * Returns a Jmol script that rotates a structure without symmetry (C1)
	 * animates the rotation around the Cn axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	private String getJmolAnimationC1(int delay) {
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

		// show inertia axes
		s.append(jmolDrawInertiaAxes());

		RectangularPrism p = new RectangularPrism(zMax-zMin, xMax-xMin, yMax-yMin);
		s.append(JmolGeometricObjects.getJmolWirePolyhedron(100, p, calcPrismTransformation()));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show Side view 1
		s.append("set echo top center;");
		s.append("echo Top view;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Back view
		s.append("set echo top center;");
		s.append("echo Side view;");
		s.append("move  0 90 0 0 0 0 0 0 2;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view

		s.append("set echo top center;");
		s.append("echo;");
		s.append("move 0 -90 0 0 0 0 0 0 2;");
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("move -90 0 0 0 0 0 0 0 2;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}
	/**
	 * Returns a Jmol script that rotates a structure with
	 * cyclic (Cn) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the Cn axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	private String getJmolAnimationC2(int delay) {
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

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo ");
		s.append("2-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
		s.append(getJmolSymmetryAxes());

		// TODO the center for Cn should be the geometric center, not the centroid (part of revereseTransformationMatrix)

		RectangularPrism p = new RectangularPrism(getZRadius()*2, getXRadius()*2, getYRadius()*2);
		s.append(JmolGeometricObjects.getJmolWirePolyhedron(100, p, calcPrismTransformation()));

		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around 2-fold axis
		s.append("move 0 180 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Side view 1
		s.append("set echo top center;");
		s.append("echo Top view ;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Back view
		s.append("set echo top center;");
		s.append("echo Side view;");
		s.append("move 0 90 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo;");
		s.append("move 0 -90 0 0 0 0 0 0 2;");
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("move -90 0 0 0 0 0 0 0 2;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
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

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo ");
		s.append(fold);
		s.append("-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
//		s.append(jmolDrawInertiaAxes());
		s.append(getJmolSymmetryAxes());

		// TODO the center for Cn should be the geometric center, not the centroid (part of revereseTransformationMatrix)
	
		Prism p = new Prism(fold);
		p.setHeight(getZRadius()*2);
		p.setInscribedRadius(getXYRadius());
		s.append(JmolGeometricObjects.getJmolWirePolyhedron(100, p, calcPrismTransformation()));
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

		Prism p = new Prism(fold);
		p.setHeight(getZRadius()*2);
		p.setInscribedRadius(getXYRadius());
		s.append(JmolGeometricObjects.getJmolWirePolyhedron(100, p, calcPrismTransformation()));
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
		Tetrahedron t = new Tetrahedron();
		t.setMidRadius(radius);
		s.append(JmolGeometricObjects.getJmolWirePolyhedron(100, t, reverseTransformationMatrix));
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
//		s.append(getJmolSymmetryAxes());
		double radius = Math.max(getZRadius(), getXYRadius());
		Octahedron o = new Octahedron();
		o.setMidRadius(radius);
		s.append(JmolGeometricObjects.getJmolWirePolyhedron(100, o, reverseTransformationMatrix));
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
//		s.append(getJmolSymmetryAxes());
		double radius = Math.max(getZRadius(), getXYRadius());
		Icosahedron p = new Icosahedron();
		p.setMidRadius(radius);
		s.append(JmolGeometricObjects.getJmolWirePolyhedron(100, p, reverseTransformationMatrix));
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
		s.append("move 100.81231696357165 0 0 0 0 0 0  0 4;");
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
		s.append("move 21 0 0 0 0 0 0 0 4;");
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
		s.append("move -121.812 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

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
