package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.quaternary.geometry.MomentsOfInertia;
import org.biojava3.structure.quaternary.geometry.SuperPosition;

public class HelixAxisAligner extends AxisAligner {
	private static final Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static final Vector3d Z_AXIS = new Vector3d(0,0,1);

	private Subunits subunits = null;
	private HelixLayers helixLayers = null;

	private Matrix4d transformationMatrix = new Matrix4d();
	private Matrix4d reverseTransformationMatrix = new Matrix4d();
	private Vector3d referenceVector = new Vector3d();
	private Vector3d principalRotationVector = new Vector3d();
	private Vector3d[] principalAxesOfInertia = null;
	private List<List<Integer>> alignedOrbits = null;

	private Vector3d minBoundary = new Vector3d();
	private Vector3d maxBoundary = new Vector3d();
	private double xyRadiusMax = Double.MIN_VALUE;

	boolean modified = true;

	public HelixAxisAligner(QuatSymmetryResults results) {
		this.subunits = results.getSubunits();
		this.helixLayers = results.getHelixLayers();
		if (subunits == null) {
			throw new IllegalArgumentException("HelixAxisTransformation: Subunits are null");
		} else if (helixLayers == null) {
			throw new IllegalArgumentException("HelixAxisTransformation: HelixLayers is null");
		} else if (subunits.getSubunitCount() == 0) {
			throw new IllegalArgumentException("HelixAxisTransformation: Subunits is empty");
		} else if (helixLayers.size() == 0) {
			throw new IllegalArgumentException("HelixAxisTransformation: HelixLayers is empty");
		}
	}

	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getTransformation()
	 */
	@Override
	public String getSymmetry() {
		run();   
		return "H";
	}
	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getTransformation()
	 */
	@Override
	public Matrix4d getTransformation() {
		run();   
		return transformationMatrix;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getRotationMatrix()
	 */
	@Override
	public Matrix3d getRotationMatrix() {
		run();
		Matrix3d m = new Matrix3d();
		transformationMatrix.getRotationScale(m);
	    return m;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getReverseTransformation()
	 */
	@Override
	public Matrix4d getReverseTransformation() {
		run();   
		return reverseTransformationMatrix;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getPrincipalRotationAxis()
	 */
	@Override
	public Vector3d getPrincipalRotationAxis() {
		run();
		return principalRotationVector;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getRotationReferenceAxis()
	 */
	@Override
	public Vector3d getRotationReferenceAxis() {
		run();
		return referenceVector;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getPrincipalAxesOfInertia()
	 */
	@Override
	public Vector3d[] getPrincipalAxesOfInertia() {
		run();
		return principalAxesOfInertia;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getDimension()
	 */
	@Override
	public Vector3d getDimension() {
		run();
		Vector3d dimension = new Vector3d();
		dimension.sub(maxBoundary, minBoundary);
		dimension.scale(0.5);
		return dimension;
	}

	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getXYRadius()
	 */
	@Override
	public double getXYRadius() {
		run();
		return xyRadiusMax;
	}

	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getGeometicCenterTransformation()
	 */
	@Override
	public Matrix4d getGeometicCenterTransformation() {
		run();
	
	    Matrix4d geometricCentered = new Matrix4d(reverseTransformationMatrix);
	    geometricCentered.setTranslation(new Vector3d(getGeometricCenter()));
	
		return geometricCentered;
	}

	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getGeometricCenter()
	 */
	@Override
	public Point3d getGeometricCenter() {
		run();
		
		Point3d geometricCenter = new Point3d();
		Vector3d translation = new Vector3d();
		reverseTransformationMatrix.get(translation);
		
		// TODO does this apply to the helic case?
		// calculate adjustment around z-axis and transform adjustment to
		//  original coordinate frame with the reverse transformation

		Vector3d corr = new Vector3d(0,0, minBoundary.z+getDimension().z);
		reverseTransformationMatrix.transform(corr);
		geometricCenter.set(corr);

		geometricCenter.add(translation);
		return geometricCenter;
	}

	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getCentroid()
	 */
	@Override
	public Point3d getCentroid() {
		return new Point3d(subunits.getCentroid());
	}

	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getSubunits()
	 */
	@Override
	public Subunits getSubunits() {
		return subunits;
	}

	public HelixLayers getHelixLayers() {
		run();
		return helixLayers;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.core.AxisAligner#getOrbits()
	 */
	@Override
	public List<List<Integer>> getOrbits() {
		run();
		return alignedOrbits;
	}

	/**
	 * @return
	 */
	
	private void run () {
		if (modified) {
			calcPrincipalRotationVector();
			calcPrincipalAxes();	
			calcReferenceVector();
			calcTransformation();
			calcReverseTransformation();
			calcBoundaries();
			calcAlignedOrbits();
			
			// orient helix along Y axis by rotating 90 degrees around X-axis
			transformationMatrix = reorientHelix(0);
			calcReverseTransformation();

			modified = false;
		}
	}

	private Matrix4d reorientHelix(int index) {
		Matrix4d matrix = new Matrix4d();
		matrix.setIdentity();
		matrix.setRotation(new AxisAngle4d(1,0,0,Math.PI/2*(index+1)));
		matrix.mul(transformationMatrix);
        return matrix;
	}
	
	/**
	 * Returns a list of orbits (an orbit is a cyclic permutation of subunit indices that are related 
	 * by a rotation around the principal rotation axis) ordered from the +z direction to the -z direction (z-depth).
	 * Within an orbit, subunit indices are ordered with a clockwise rotation around the z-axis.
	 * @return list of orbits ordered by z-depth
	 */
	private void calcAlignedOrbits() {
		Map<Double, List<Integer>> depthMap = new TreeMap<Double, List<Integer>>();
		double[] depth = getSubunitZDepth();
		alignedOrbits = calcOrbits();

		// calculate the mean depth of orbit along z-axis
		for (List<Integer> orbit: alignedOrbits) {
			// calculate the mean depth along z-axis for each orbit
			double meanDepth = 0;
			for (int subunit: orbit) {
				meanDepth += depth[subunit];
			}
			meanDepth /= orbit.size();

			if (depthMap.get(meanDepth) != null) {
				// System.out.println("Conflict in depthMap");
				meanDepth += 0.01;
			}
//			System.out.println("calcAlignedOrbits: " + meanDepth + " orbit: " + orbit);
			depthMap.put(meanDepth, orbit);
		}

		// now fill orbits back into list ordered by depth
		alignedOrbits.clear();
		for (List<Integer> orbit: depthMap.values()) {
			// order subunit in a clockwise rotation around the z-axis 
			/// starting at the 12 O-clock position (+y position)
			// TODO how should this be aligned??
	//		alignWithReferenceAxis(orbit);
			alignedOrbits.add(orbit);
		}
	}
	
	

	

	private void calcTransformation() {
		calcTransformationBySymmetryAxes();
		// make sure this value is zero. On Linux this value is 0.
		transformationMatrix.setElement(3, 3, 1.0);
	}

	private void calcReverseTransformation() {
		reverseTransformationMatrix.invert(transformationMatrix);
	}

	private void calcTransformationBySymmetryAxes() {
		Vector3d[] axisVectors = new Vector3d[2];
		axisVectors[0] = new Vector3d(principalRotationVector);
		axisVectors[1] = new Vector3d(referenceVector);

		//  y,z axis centered at the centroid of the subunits
		Vector3d[] referenceVectors = new Vector3d[2];
		referenceVectors[0] = new Vector3d(Z_AXIS);
		referenceVectors[1] = new Vector3d(Y_AXIS);

		transformationMatrix = alignAxes(axisVectors, referenceVectors);

		// combine with translation
		Matrix4d combined = new Matrix4d();
		combined.setIdentity();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		combined.setTranslation(trans);
		transformationMatrix.mul(combined);
		
		// for helical geometry, set a canonical view for the Z direction
		calcZDirection();
	}

	/**
	 * Returns a transformation matrix that rotates refPoints to match
	 * coordPoints
	 * @param refPoints the points to be aligned
	 * @param referenceVectors
	 * @return
	 */
	private Matrix4d alignAxes(Vector3d[] axisVectors, Vector3d[] referenceVectors) {
		Matrix4d m1 = new Matrix4d();
		AxisAngle4d a = new AxisAngle4d();
		Vector3d axis = new Vector3d();
		
		// calculate rotation matrix to rotate refPoints[0] into coordPoints[0]
		Vector3d v1 = new Vector3d(axisVectors[0]);
		Vector3d v2 = new Vector3d(referenceVectors[0]);
		double dot = v1.dot(v2);
		if (Math.abs(dot) < 0.999) {
			axis.cross(v1,v2);
			axis.normalize();
			a.set(axis, v1.angle(v2));
			m1.set(a);
			// make sure matrix element m33 is 1.0. It's 0 on Linux.
			m1.setElement(3,  3, 1.0);
		} else if (dot > 0) {
			// parallel axis, nothing to do -> identity matrix
			m1.setIdentity();
		} else if (dot < 0) {
			// anti-parallel axis, flip around x-axis
			m1.set(flipX());
		}
		
		// apply transformation matrix to all refPoints
		m1.transform(axisVectors[0]);
		m1.transform(axisVectors[1]);
		
		// calculate rotation matrix to rotate refPoints[1] into coordPoints[1]
		v1 = new Vector3d(axisVectors[1]);
		v2 = new Vector3d(referenceVectors[1]);
		Matrix4d m2 = new Matrix4d();
		dot = v1.dot(v2);
		if (Math.abs(dot) < 0.999) {
			axis.cross(v1,v2);
			axis.normalize();
			a.set(axis, v1.angle(v2));
			m2.set(a);
			// make sure matrix element m33 is 1.0. It's 0 on Linux.
			m2.setElement(3,  3, 1.0);
		} else if (dot > 0) {
			// parallel axis, nothing to do -> identity matrix
			m2.setIdentity();
		} else if (dot < 0) {
			// anti-parallel axis, flip around z-axis
			m2.set(flipZ());
		}
		
		// apply transformation matrix to all refPoints
		m2.transform(axisVectors[0]);
		m2.transform(axisVectors[1]);
		
		// combine the two rotation matrices
		m2.mul(m1);

		// the RMSD should be close to zero
		Point3d[] axes = new Point3d[2];
		axes[0] = new Point3d(axisVectors[0]);
		axes[1] = new Point3d(axisVectors[1]);
		Point3d[] ref = new Point3d[2];
		ref[0] = new Point3d(referenceVectors[0]);
		ref[1] = new Point3d(referenceVectors[1]);
		if (SuperPosition.rmsd(axes, ref) > 0.1) {
			System.out.println("Warning: AxisTransformation: axes alignment is off. RMSD: " + SuperPosition.rmsd(axes, ref));
		}
		
		return m2;
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

	/**
	 * Calculates the min and max boundaries of the structure after it has been
	 * transformed into its canonical orientation.
	 */
	private void calcBoundaries() {	
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
	}
	
	/*
	 * Modifies the rotation part of the transformation axis for
	 * a Cn symmetric complex, so that the narrower end faces the
	 * viewer, and the wider end faces away from the viewer. Example: 3LSV
	 */
	private void calcZDirection() {
		calcBoundaries();
	
		// if the longer part of the structure faces towards the back (-z direction),
		// rotate around y-axis so the longer part faces the viewer (+z direction)
		if (Math.abs(minBoundary.z) > Math.abs(maxBoundary.z)) {
			Matrix4d rot = flipY();
			rot.mul(transformationMatrix);
			transformationMatrix.set(rot);
		}
	}

	private double[] getSubunitZDepth() {	
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

		List<List<Integer>> orbits = new ArrayList<List<Integer>>();
		for (int i = 0; i < n; i++) {
			orbits.add(Collections.singletonList(i));
		}
		
		return orbits;
	}
		
	/**
	 * Returns a vector along the principal rotation axis for the
	 * alignment of structures along the z-axis
	 * @return principal rotation vector
	 */
	private void calcPrincipalRotationVector() {
		AxisAngle4d axisAngle = helixLayers.getByLowestAngle().getAxisAngle();
		principalRotationVector = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
	}

	/**
	 * Returns a vector perpendicular to the principal rotation vector
	 * for the alignment of structures in the xy-plane
	 * @return reference vector
	 */
	private void calcReferenceVector() {
		referenceVector = getReferenceAxisCylic();
		
		if (referenceVector == null) {
			System.err.println("Warning: no reference vector found. Using y-axis.");
			referenceVector = new Vector3d(Y_AXIS);
		}
		// make sure reference vector is perpendicular principal roation vector
		referenceVector = orthogonalize(principalRotationVector, referenceVector);
	}
	
	private Vector3d orthogonalize(Vector3d vector1, Vector3d vector2) {
		double dot = vector1.dot(vector2);
		Vector3d ref = new Vector3d(vector2);
//		System.out.println("p.r: " + dot);
//		System.out.println("Orig refVector: " + referenceVector);
		if (dot < 0) {
			vector2.negate();
		}
		vector2.cross(vector1, vector2);
//		System.out.println("Intermed. refVector: " + vector2);
		vector2.normalize();
//		referenceVector.cross(referenceVector, principalRotationVector); 
		vector2.cross(vector1, vector2); 
		vector2.normalize();	
		if (ref.dot(vector2) < 0) {
			vector2.negate();
		}
//		System.out.println("Mod. refVector: " + vector2);
		return vector2;
	}
	/**
	 * Returns the default reference vector for the alignment of Cn structures
	 * @return
	 */	
	private Vector3d getReferenceAxisCylic() {		
		// get principal axis vector that is perpendicular to the principal 
		// rotation vector
		Vector3d vmin = null;
		double dotMin = 1.0;
		for (Vector3d v: principalAxesOfInertia) {
			if (Math.abs(principalRotationVector.dot(v)) < dotMin) {
				dotMin = Math.abs(principalRotationVector.dot(v));
				vmin = new Vector3d(v);
			}
		}
		if (principalRotationVector.dot(vmin) < 0) {
			vmin.negate();
		}
		
		return vmin;
	}

	private static Matrix4d flipX() {
		Matrix4d rot = new Matrix4d();
		rot.m00 = 1;
		rot.m11 = -1;
		rot.m22 = -1;
		rot.m33 = 1;
		return rot;
	}
	
	private static Matrix4d flipY() {
		Matrix4d rot = new Matrix4d();
		rot.m00 = -1;
		rot.m11 = 1;
		rot.m22 = -1;
		rot.m33 = 1;
		return rot;
	}
	
	private static Matrix4d flipZ() {
		Matrix4d rot = new Matrix4d();
		rot.m00 = -1;
		rot.m11 = -1;
		rot.m22 = 1;
		rot.m33 = 1;
		return rot;
	}
	
}
