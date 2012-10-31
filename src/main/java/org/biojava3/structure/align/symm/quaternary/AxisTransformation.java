package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

public class AxisTransformation {
	private static Vector3d X_AXIS = new Vector3d(1,0,0);
	private static Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static Vector3d Z_AXIS = new Vector3d(0,0,1);

	private Subunits subunits = null;
	private RotationGroup rotationGroup = null;

	private Matrix4d transformationMatrix = new Matrix4d();
	private Matrix4d reverseTransformationMatrix = new Matrix4d();
	private Vector3d principalRotationAxis = new Vector3d();
	private Vector3d referenceAxis = new Vector3d();
	private Vector3d[] principalAxesOfInertia = null;

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
		run();   
		return transformationMatrix;
	}
	
	public Matrix4d getReverseTransformation() {
		run();   
		return reverseTransformationMatrix;
	}
	
	public Vector3d getPrincipalRotationAxis() {
		return principalRotationAxis;
	}
	
	public Vector3d getRotationReferenceAxis() {
		return referenceAxis;
	}
	
	public Vector3d[] getPrincipalAxesOfInertia() {
		return principalAxesOfInertia;
	}
	
	/**
	 * Returns the x-radius
	 * @return double radius in z-direction
	 */
	public double getXRadius() {
		return 0.5 * (xMax - xMin); // half of dimension along z-axis (principal rotation axis)
	}
	/**
	 * Returns the y-radius
	 * @return double radius in z-direction
	 */
	public double getYRadius() {
		return 0.5 * (yMax - yMin); // half of dimension along z-axis (principal rotation axis)
	}
	
	/**
	 * Returns the z-radius
	 * @return double radius in z-direction
	 */
	public double getZRadius() {
		return 0.5 * (zMax - zMin); // half of dimension along z-axis (principal rotation axis)
	}
	
	/**
	 * Returns the radius for drawing the minor rotation axis in the xy-plane
	 * @return double radius in xy-plane
	 */
	public double getXYRadius() {
		return xyRadiusMax;
	}


	private void run () {
		if (modified) {
			calcPrincipalRotationAxis();
			calcPrincipalAxes();
			calcReferenceAxis();
			calcTransformation();
			calcBoundaries();
			modified = false;
		}
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
	//	principalAxis = getPrincipalRotationAxis();
		calcReferenceAxis();

		Point3d[] refPoints = new Point3d[2];
		refPoints[0] = new Point3d(principalRotationAxis);
		refPoints[1] = new Point3d(referenceAxis);


		//  y,z axis centered at the centroid of the subunits
		Point3d[] coordPoints = new Point3d[2];
		coordPoints[0] = new Point3d(Z_AXIS);
		coordPoints[1] = new Point3d(Y_AXIS);

		Matrix4d matrix = alignAxes(refPoints, coordPoints);

		calcReverseMatrix(matrix);

		return matrix;
	}

	private void calcReferenceAxis() {
		referenceAxis = getMinorRotationAxis();
		if (referenceAxis == null) {
			referenceAxis = getSubunitReferenceAxisNew();
		}
		// make sure reference axis is perpendicular
		referenceAxis.cross(principalRotationAxis, referenceAxis);
		referenceAxis.cross(referenceAxis, principalRotationAxis);
	}
	
	private Matrix4d getTransformationByC2SymmetryAxes() {
	calcReferenceAxis();
		
		Point3d[] refPoints = new Point3d[2];
		refPoints[0] = new Point3d(principalRotationAxis);
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
		
		Vector3d vmin = null;
		double dotMin = 1.0;
		for (Vector3d v: principalAxesOfInertia) {
			if (Math.abs(principalRotationAxis.dot(v)) < dotMin) {
				dotMin = Math.abs(principalRotationAxis.dot(v));
				vmin = new Vector3d(v);
			}
		}
		return vmin;
	}


	/**
	 * Returns a transformation matrix for prisms to draw for Cn structures.
	 * The center in this matrix is the geometric center, rather then the centroid
	 * In Cn structures those are usually not the same.
	 * @return
	 */
	public Matrix4d calcPrismTransformation() {
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

	public Point3d calcGeometricCenter() {
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
	private void calcPrincipalRotationAxis() {
		Rotation rotation = rotationGroup.getRotation(0); // the rotation around the principal axis is the first rotation
		AxisAngle4d axisAngle = rotation.getAxisAngle();
		principalRotationAxis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		principalRotationAxis.normalize();
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
				double d = v.dot(principalRotationAxis);
				if (rotationGroup.getRotation(i).getFold() == 5) {
					// the dot product of 0.447.. is between to adjacent 5-fold axis
					if (d > 0.447 && d < 0.448) {
						return v;
					}
				}
		}
		return null;
	}

	private void calcPrincipalAxes() {
		MomentsOfInertia moi = new MomentsOfInertia();

		for (Point3d[] list: subunits.getTraces()) {
			for (Point3d p: list) {
				moi.addPoint(p, 1.0);
			}
		}
		principalAxesOfInertia = moi.getPrincipalAxes();
	}

	private Matrix4d getTransformationByInertiaAxes() {
		Point3d[] refPoints = new Point3d[2];
		refPoints[0] = new Point3d(principalAxesOfInertia[0]);
		refPoints[1] = new Point3d(principalAxesOfInertia[1]);

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

}
