package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.quaternary.geometry.MomentsOfInertia;
import org.biojava3.structure.quaternary.geometry.SuperPosition;

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

	private Vector3d minBoundary = new Vector3d();
	private Vector3d maxBoundary = new Vector3d();
	private double xyRadiusMax = Double.MIN_VALUE;

	boolean modified = true;

	public AxisTransformation(Subunits subunits, RotationGroup rotationGroup) {
		this.subunits = subunits;
		this.rotationGroup = rotationGroup;
		if (subunits == null) {
			throw new IllegalArgumentException("AxisTransformation: Subunits are null");
		} else if (rotationGroup == null) {
			throw new IllegalArgumentException("AxisTransformation: RotationGroup is null");
		} else if (subunits.getSubunitCount() == 0) {
			throw new IllegalArgumentException("AxisTransformation: Subunits is empty");
		} else if (rotationGroup.getOrder() == 0) {
			throw new IllegalArgumentException("AxisTransformation: RotationGroup is empty");
		}
	}

	public Matrix4d getTransformation() {
		run();   
		return transformationMatrix;
	}
	
	public Matrix3d getRotationMatrix() {
		run();
		Matrix3d m = new Matrix3d();
		transformationMatrix.getRotationScale(m);
	    return m;
	}
	
	public Matrix4d getReverseTransformation() {
		run();   
		return reverseTransformationMatrix;
	}
	
	public Vector3d getPrincipalRotationAxis() {
		run();
		return principalRotationAxis;
	}
	
	public Vector3d getRotationReferenceAxis() {
		run();
		return referenceAxis;
	}
	
	public Vector3d[] getPrincipalAxesOfInertia() {
		run();
		return principalAxesOfInertia;
	}
	
	public Vector3d getDimension() {
		run();
		Vector3d dimension = new Vector3d();
		dimension.sub(maxBoundary, minBoundary);
		dimension.scale(0.5);
		return dimension;
	}

	/**
	 * Returns the radius for drawing the minor rotation axis in the xy-plane
	 * @return double radius in xy-plane
	 */
	public double getXYRadius() {
		run();
		return xyRadiusMax;
	}

	public Subunits getSubunits() {
		return subunits;
	}

	public RotationGroup getRotationGroup() {
		return rotationGroup;
	}
	
	private void run () {
		if (modified) {
			calcPrincipalRotationAxis();
			calcPrincipalAxes();
			
			// initial alignment with draft ref. axis
			boolean refine = false;
			calcReferenceAxis(refine);
			calcTransformation();
			calcReverseTransformation();
			calcBoundaries();
			
			// refine ref. axis for cyclic systems
			if (rotationGroup.getPointGroup().startsWith("C") &&
					!rotationGroup.getPointGroup().startsWith("C2")) {
				refine = true;
				calcReferenceAxis(refine);
				calcTransformation();
				calcReverseTransformation();
				calcBoundaries();
			}
			modified = false;
		}
	}

	private void calcTransformation() {
		if (subunits == null || subunits.getSubunitCount() == 0) {
			// if the structure doesn't contain protein subunits, i.e.,
			// a DNA structure, return the identity matrix.
			transformationMatrix.setIdentity();
		} else if (rotationGroup.getPointGroup().equals("C1")) {
			transformationMatrix = getTransformationByInertiaAxes();
//		} else if (rotationGroup.getPointGroup().equals("C2")) {
//			transformationMatrix = getTransformationByC2SymmetryAxes();;
		} else {
			transformationMatrix = getTransformationBySymmetryAxes();
			// for cyclic geometry, set a canonical view for the Z direction
			if (rotationGroup.getPointGroup().startsWith("C")) {
				setZDirection(transformationMatrix);
			}
		}
	}

	private Matrix4d getTransformationBySymmetryAxes() {
		Point3d[] refPoints = new Point3d[2];
		refPoints[0] = new Point3d(principalRotationAxis);
		refPoints[1] = new Point3d(referenceAxis);
//		System.out.println("getTransfromationBySymmetryAxes: refAxis:" + referenceAxis);


		//  y,z axis centered at the centroid of the subunits
		Point3d[] coordPoints = new Point3d[2];
		coordPoints[0] = new Point3d(Z_AXIS);
		coordPoints[1] = new Point3d(Y_AXIS);

		Matrix4d matrix = alignAxes(refPoints, coordPoints);

		// combine with translation
		Matrix4d combined = new Matrix4d();
		combined.setIdentity();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		combined.setTranslation(trans);
		matrix.mul(combined);
	    return matrix;
	}

	private void calcReferenceAxis(boolean refine) {
		referenceAxis = getReferenceAxis(refine);
//		if (referenceAxis == null) {
//			System.out.println("Ref. axis null");
//			if (refine) {
//				referenceAxis = getSubunitReferenceAxisByOrbit();
//			} else {
//				referenceAxis = getSubunitReferenceAxis();
//			}
//		}
		// make sure reference axis is perpendicular
		referenceAxis.cross(principalRotationAxis, referenceAxis);
		referenceAxis.cross(referenceAxis, principalRotationAxis); 
	}
	
	private Matrix4d getTransformationByC2SymmetryAxes() {		
		Point3d[] refPoints = new Point3d[2];
		refPoints[0] = new Point3d(principalRotationAxis);
		refPoints[1] = new Point3d(referenceAxis);
		
		Point3d[] coordPoints = new Point3d[2];
		coordPoints[0] = new Point3d(Y_AXIS);
		coordPoints[1] = new Point3d(X_AXIS);
 
		Matrix4d matrix = alignAxes(refPoints, coordPoints);  
		
		// combine with translation
		Matrix4d translation = new Matrix4d();
		translation.setIdentity();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		translation.setTranslation(trans);
		matrix.mul(translation);
		return matrix;
	}
	
	private Matrix4d alignAxes(Point3d[] refPoints, Point3d[] coordPoints) {
		Vector3d v1 = new Vector3d(refPoints[0]);
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
		
		v1 = new Vector3d(refPoints[1]);
		v2 = new Vector3d(coordPoints[1]);
		axis.cross(v1,v2);
		a.set(axis, v1.angle(v2));
		Matrix4d m2 = new Matrix4d();
		m2.set(a);
		for (Point3d v: refPoints) {
			m2.transform(v);
		}

		if (SuperPosition.rmsd(refPoints, coordPoints) > 0.01) {
			System.out.println("Warning: AxisTransformation: axes alignment is off. RMSD: " + SuperPosition.rmsd(refPoints, coordPoints));
		}
		m2.mul(m1);
		return m2;
	}

	private void calcReverseTransformation() {
		reverseTransformationMatrix.invert(transformationMatrix);
	}


	/**
	 * Returns a transformation matrix transform polyhedra for Cn structures.
	 * The center in this matrix is the geometric center, rather then the centroid.
	 * In Cn structures those are usually not the same.
	 * @return
	 */
	public Matrix4d getGeometicCenterTransformation() {
		run();

	    Matrix4d geometricCentered = new Matrix4d(reverseTransformationMatrix);
	    geometricCentered.setTranslation(new Vector3d(getGeometricCenter()));

		return geometricCentered;
	}

	/**
	 * Returns the geometric center of polyhedron. In the case of the Cn 
	 * point group, the centroid and geometric center are usually not
	 * identical.
	 * @return
	 */
	public Point3d getGeometricCenter() {
		run();
		
		Point3d geometricCenter = new Point3d();
		Vector3d translation = new Vector3d();
		reverseTransformationMatrix.get(translation);
		
		// calculate adjustment around z-axis and transform adjustment to
		//  original coordinate frame with the reverse transformation
		if (rotationGroup.getPointGroup().startsWith("C")) {
			Vector3d corr = new Vector3d(0,0, minBoundary.z+getDimension().z);
			reverseTransformationMatrix.transform(corr);
			geometricCenter.set(corr);
		}
		
		geometricCenter.add(translation);
		return geometricCenter;
	}
	
	public Point3d getCentroid() {
		return new Point3d(subunits.getCentroid());
	}
	
	/**
	 * Calculates the min and max boundaries of the structure after it has been
	 * transformed into its canonical orientation.
	 */
	private void calcBoundaries() {
		if (subunits == null || subunits.getSubunitCount() == 0) {
			return;
		}
		
		minBoundary.x = Double.MAX_VALUE;
		maxBoundary.x = Double.MIN_VALUE;
		minBoundary.y = Double.MAX_VALUE;
		maxBoundary.x = Double.MIN_VALUE;
		minBoundary.z = Double.MAX_VALUE;
		maxBoundary.z = Double.MIN_VALUE;

		Point3d probe = new Point3d();

		for (Point3d[] list: subunits.getTraces()) {
			for (Point3d p: list) {
				probe.set(p);
				transformationMatrix.transform(probe);

				minBoundary.x = Math.min(minBoundary.x, probe.x);
				maxBoundary.x = Math.max(maxBoundary.x, probe.x);
				minBoundary.y = Math.min(minBoundary.y, probe.y);
				maxBoundary.y = Math.max(maxBoundary.y, probe.y);
				minBoundary.z = Math.min(minBoundary.z, probe.z);
				maxBoundary.z = Math.max(maxBoundary.z, probe.z);
				xyRadiusMax = Math.max(xyRadiusMax, Math.sqrt(probe.x*probe.x + probe.y * probe.y));
			}
		}
//		System.out.println("Min: " + minBoundary);
//		System.out.println("Max: " + maxBoundary);
	}
	
	/**
	 * 
	 */
	public List<List<Integer>> getOrbits() {
		Map<Double, List<Integer>> depthMap = new TreeMap<Double, List<Integer>>();
		double[] depth = getOrbitDepth();
		List<List<Integer>> orbits = calcOrbits();

		// calculate the mean depth of orbit along z-axis
		for (List<Integer> orbit: orbits) {
			// calculate the mean depth along z-axis for each orbit
			double meanDepth = 0;
			for (int subunit: orbit) {
				meanDepth += depth[subunit];
			}
		    meanDepth /= orbit.size();
		    
		    if (depthMap.get(meanDepth) != null) {
//		    	System.out.println("Conflict in depthMap");
		    	meanDepth += 0.01;
		    }
		    depthMap.put(meanDepth, orbit);
		}
		
		// now fill orbits back into list ordered by depth
		orbits.clear();
		for (List<Integer> orbit: depthMap.values()) {
			orbits.add(orbit);
		}
		return orbits;
	}
	
	private double[] getOrbitDepth() {	
		int n = subunits.getSubunitCount();
        double[] depth = new double[n];        
		Point3d probe = new Point3d();

		// transform subunit centers into z-aligned position and calculate
		// z-coordinates (depth) along the z-axis.
		for (int i = 0; i < n; i++) {
			Point3d p= subunits.getCenters().get(i);
			probe.set(p);
			transformationMatrix.transform(probe);
			depth[i] = probe.z;
		}
		return depth;
	}
	
	/**
	 * Returns a list of list of subunit ids that form an "orbit", i.e. they
	 * are transformed into each other during a rotation around the principal symmetry axis (z-axis)
	 * @return
	 */
	private List<List<Integer>> calcOrbits() {
		int n = subunits.getSubunitCount();
//		for (int i = 0; i < n; i++) {
//			System.out.println("ChainId: " + subunits.getChainIds().get(i));
//		}
		int fold = rotationGroup.getRotation(0).getFold();

		List<List<Integer>> orbits = new ArrayList<List<Integer>>();
		Set<Integer> used = new HashSet<Integer>();
//		for (int j = 0; j < fold; j++) {
//			System.out.println("Rotation angle: " +rotationGroup.getRotation(j).getAxisAngle());
//			System.out.println(rotationGroup.getRotation(j).getPermutation());
//		}
		
		List<Integer> inOrder = rotationGroup.getRotation(0).getPermutation();
		
		// for simple Cn group, order the orbits in rotation order for coloring
		if (rotationGroup.getOrder() > 1  && n == fold && rotationGroup.getPointGroup().startsWith("C")) {
		    inOrder = deconvolute();
		}
			
		for (int i = 0; i < n; i++) {
			if (! used.contains(i)) {
				// determine the equivalent subunits
				List<Integer> orbitElements = new ArrayList<Integer>();
				for (int j = 0; j < fold; j++) {
					List<Integer> permutation = rotationGroup.getRotation(j).getPermutation();
					orbitElements.add(permutation.get(i));
					used.add(permutation.get(i));
				}

//				System.out.println("OrbitElements: " + orbitElements);

				// order subunits in rotation order
				List<Integer> orbit = new ArrayList<Integer>(orbitElements.size());
				for (Integer p: inOrder) {
					if (orbitElements.contains(p)) {
						orbit.add(p);
					}
				}
				orbits.add(orbit);
			}
		}
//		System.out.println("Orbits: ");
//		for (List<Integer> orbit: orbits) {
//			System.out.println(orbit);
//		}
		return orbits;
	}
	
	private List<Integer> deconvolute() {
		// get first rotation in rotation group (by definition the first rotation with smallest angle)
		List<Integer> permutation = rotationGroup.getRotation(1).getPermutation();

		List<Integer> inRotationOrder = new ArrayList<Integer>(permutation.size());
		inRotationOrder.add(0);
		for (int i = 0; i < permutation.size()-1; i++) {
			int next = permutation.get(inRotationOrder.get(i));
			if (next == 0) {
				next = i+1;
			}
		    inRotationOrder.add(next);
			
//			System.out.println("inrotationorder: " + inRotationOrder);
		}
		return inRotationOrder;
	}
	
	/*
	 * Modifies the rotation part of the transformation axis for
	 * a Cn symmetric complex, so that the narrower end faces the
	 * viewer, and the wider end faces away from the viewer. Example: 3LSV
	 */
	private Matrix4d setZDirection(Matrix4d matrix) {
		calcBoundaries();

		// if the longer part of the structure faces towards the back (-z direction),
		// rotate around y-axis so the longer part faces the viewer (+z direction)
		if (Math.abs(minBoundary.z) > Math.abs(maxBoundary.z)) {
			Matrix4d rot = new Matrix4d();
			rot.rotY(-Math.PI);
			rot.mul(matrix);
			matrix.set(rot);
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
	 * Returns a normalized vector that represents the minor rotation axis, except 
	 * for Cn, this represents an axis orthogonal to the principal axis
	 * @return minor rotation axis, return null if minor rotation axis does not exists
	 */
	private Vector3d getReferenceAxis(boolean refine) {
		if (rotationGroup.getPointGroup().startsWith("C")) {
			if (refine) {
				return getSubunitReferenceAxisByOrbit();
			} else {
				return getSubunitReferenceAxis();
			}
		} else if (rotationGroup.getPointGroup().startsWith("D")) {
			return getReferenceAxisDihedral();
		}else if (rotationGroup.getPointGroup().equals("T")) {
			return getReferenceAxisTetrahedral();
		} else if (rotationGroup.getPointGroup().equals("O")) {
			return getReferenceAxisOctahedral();
		} else if (rotationGroup.getPointGroup().equals("I")) {
			return getReferenceAxisIcosahedral();		
		} else {
			return null; 
		}
	}
	
	/**
	 * Returns vector from largest subunit to centroid of complex
	 * @return
	 */	
	private Vector3d getSubunitReferenceAxis() {
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

	private Vector3d getSubunitReferenceAxisByOrbit() {
		if (rotationGroup.getPointGroup().equals("C2")) {
			return referenceAxis;
		} 

		Point3d probe = new Point3d();
		double xyMax = Double.MIN_VALUE;
		
		List<Point3d> centers = subunits.getCenters();
		int refIndex = 0;
		for (int i = 0; i < centers.size(); i++) {
			probe.set(centers.get(i));
			transformationMatrix.transform(probe);
			double xy = probe.x*probe.x + probe.y*probe.y;
			if (xy > xyMax) {
				xyMax = xy;
				refIndex = i;
			}
		}
		Vector3d refAxis = new Vector3d();
		refAxis.sub(subunits.getCentroid(), centers.get(refIndex));
		refAxis.normalize();
		return refAxis;
	}


	/**
	 * 
	 */
	private Vector3d getReferenceAxisDihedral() {
		System.out.println("getReferenceAxisDihedral");
		int maxFold = rotationGroup.getRotation(0).getFold();
		// one exception: D2
		if (maxFold == 2) {
			maxFold = 3;
		}
    	// TODO how about D2, where minor axis = 2 = principal axis??
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
			if (rotationGroup.getRotation(i).getDirection() == 1 && rotationGroup.getRotation(i).getFold() < maxFold) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				v.normalize();
				// for consistency with the other getReferenceAxis... method, use positive dot product
				// to provide uniform directionality
				if (v.dot(principalRotationAxis) < 0) {
					continue;
				}
//				System.out.println("Minor rotation axis: " + rotationGroup.getRotation(i).getFold());
				return v;
			}
		}
		return null;
	}
	
	private Vector3d getReferenceAxisTetrahedral() {
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				double d = v.dot(principalRotationAxis);
				if (rotationGroup.getRotation(i).getFold() == 3) {
					// the dot product 0 is between to adjacent 3-fold axes
					if (d > 0.3 && d < 0.9) {
						return v;
					}
				}
		}
		return null;
	}
	
	private Vector3d getReferenceAxisOctahedral() {
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				double d = v.dot(principalRotationAxis);
				if (rotationGroup.getRotation(i).getFold() == 4) {
					// the dot product 0 is between to adjacent 4-fold axes
					if (d >= 0 && d < 0.1 ) {
						return v;
					}
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
					// the dot product of 0.447.. is between to adjacent 5-fold axes
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
			System.out.println("Warning: aligment with coordinate system is off. RMSD: " + SuperPosition.rmsd(refPoints, coordPoints));
		}
		
		// combine with translation
		Matrix4d translation = new Matrix4d();
		translation.setIdentity();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		translation.setTranslation(trans);
		matrix.mul(translation);
		
		return matrix;
	}

}
