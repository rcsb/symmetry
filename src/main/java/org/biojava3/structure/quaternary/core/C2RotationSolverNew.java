
package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.quaternary.geometry.CircleSampler;
import org.biojava3.structure.quaternary.geometry.MomentsOfInertia;
import org.biojava3.structure.quaternary.geometry.SuperPosition;


/**
 *
 * @author Peter
 */
public class C2RotationSolverNew implements QuatSymmetrySolver {
	private static final Vector3d X_AXIS = new Vector3d(1,0,0);
	private static final Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static final Vector3d Z_AXIS = new Vector3d(0,0,1);
	
    private Subunits subunits = null;
    private QuatSymmetryParameters parameters = null;
//    private double rmsdThreshold = 0.0f;
    private double distanceThreshold = 0.0f;
    private Vector3d centroid = new Vector3d();
    private Matrix4d centroidInverse = new Matrix4d();
    private Point3d[] originalCoords = null;
    private Point3d[] transformedCoords = null;
    private Matrix4d referenceTransformation = null;
    private double bestRmsd = Double.MAX_VALUE;

    private RotationGroup rotations = new RotationGroup();
    private QuatSuperpositionScorer scorer = null;

    public C2RotationSolverNew(Subunits subunits, QuatSymmetryParameters parameters) {
        if (subunits.getSubunitCount() != 2) {
    		throw new IllegalArgumentException("C2RotationSolver can only be applied to cases with 2 centers");
    	}
        this.subunits = subunits;
        this.parameters = parameters;
    }
   
    public RotationGroup getSymmetryOperations() {
        if (rotations.getOrder() == 0) {
            solve();
        }
        return rotations;
    }
    
    private void solve() {
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		List<Point3d[]> traces = subunits.getTraces();

		Point3d[] x = SuperPosition.clonePoint3dArray(traces.get(0));
		SuperPosition.center(x);
		Point3d[] y = SuperPosition.clonePoint3dArray(traces.get(1));
		SuperPosition.center(y);

		AxisAngle4d axisAngle = new AxisAngle4d();
		Matrix4d transformation = SuperPosition.superposeAtOrigin(x, y, axisAngle);
		System.out.println("Transformation: " + transformation);
		double caRmsd = SuperPosition.rmsd(x,  y);
		System.out.println("Rmsd: " + caRmsd);
		
		// combine with translation
		Matrix4d combined = new Matrix4d();
		combined.setIdentity();
		combined.setTranslation(trans);
		transformation.mul(combined);

		addEOperation();
		List<Integer> permutation = new ArrayList<Integer>();
		permutation.add(new Integer(1));
		permutation.add(new Integer(0));
		Rotation symmetryOperation = createSymmetryOperation(permutation, transformation, axisAngle, 0.0, caRmsd, 2);
		rotations.addRotation(symmetryOperation);
    }

    private void solve2() {
    	
    	MomentsOfInertia moi = subunits.getMomentsOfInertia();
    	Vector3d[] axes = moi.getPrincipalAxes();
    	Vector3d referenceAxis = new Vector3d();
    	referenceAxis.sub(subunits.getCenters().get(0), subunits.getCenters().get(1));
    	referenceAxis.normalize();
    	
    	System.out.println("referenceAxis: " + referenceAxis);
    	Vector3d[] axisVectors = new Vector3d[2];
    	axisVectors[0] = referenceAxis;

    	for (Vector3d v: axes) {
    		System.out.println("axisVectors: " + v);
    		if (Math.abs(referenceAxis.dot(v)) < 0.5) {
    			axisVectors[1] = orthogonalize(referenceAxis, v);
    			System.out.println("orthogonal axisVectors: " + axisVectors + " dot: " + axisVectors[0].dot(referenceAxis));
    			break;
    		}
    	}
    //	axisVectors[1] = new Vector3d();
    //	axisVectors[1].cross(referenceAxis, axisVectors[0]);

		//  y,z axis centered at the centroid of the subunits
		Vector3d[] referenceVectors = new Vector3d[2];
		referenceVectors[0] = new Vector3d(X_AXIS);
		referenceVectors[1] = new Vector3d(Z_AXIS);

		Matrix4d transformationMatrix = alignAxes(axisVectors, referenceVectors);
		transformationMatrix.transform(axisVectors[0]);
		transformationMatrix.transform(axisVectors[1]);
		System.out.println("Transformed axisVectors Z: " + axisVectors[0]);
		System.out.println("Transformed axisVectors Y: " + axisVectors[1]);
		
		

		// combine with translation
		Matrix4d combined = new Matrix4d();
		combined.setIdentity();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		combined.setTranslation(trans);
		transformationMatrix.mul(combined);
		
		double caRmsd = calcTraceRmsd(transformationMatrix);
		System.out.println("RMSD: " + caRmsd);
		addEOperation();
		List<Integer> permutation = new ArrayList<Integer>();
		permutation.add(new Integer(1));
		permutation.add(new Integer(0));
		AxisAngle4d axisAngle = new AxisAngle4d();
		axisAngle.set(transformationMatrix);
		Rotation symmetryOperation = createSymmetryOperation(permutation, transformationMatrix, axisAngle, 0.0, caRmsd, 2);
		rotations.addRotation(symmetryOperation);
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
    
    private double calcTraceRmsd(Matrix4d transformationMatrix) {
    	List<Point3d[]> traces = subunits.getTraces();
    	
    	Point3d[] trace0 = SuperPosition.clonePoint3dArray(traces.get(0));
    	SuperPosition.transform(transformationMatrix, trace0);
    	Point3d[] trace1 = SuperPosition.clonePoint3dArray(traces.get(1));
     	SuperPosition.transform(transformationMatrix, trace1);
        SuperPosition.transform(flipZ(), trace1);
    	return SuperPosition.rmsd(trace0,  trace1);
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
		if (SuperPosition.rmsd(axes, ref) > 0.01) {
			System.out.println("Warning: AxisTransformation: axes alignment is off. RMSD: " + SuperPosition.rmsd(axes, ref));
		}
		
		return m2;
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
    private void solve1() {  	
        initialize();
        
        // add the unit operation
        addEOperation();
        
        AxisAngle4d sphereAngle = new AxisAngle4d();
        Matrix4d transformation = new Matrix4d();

        int n = subunits.getSubunitCount();

        for (int i = 0; i < CircleSampler.getSphereCount(); i++) {
        	CircleSampler.getAxisAngle(i, sphereAngle);
        	transformRotationAxis(sphereAngle);    	
        	sphereAngle.angle = Math.PI;
        	transformation.set(sphereAngle);
        	// make sure matrix element m33 is 1. It's zero on Linux
        	transformation.setElement(3, 3,  1);
        	for (int j = 0; j < n; j++) {
        		transformedCoords[j].set(originalCoords[j]);
        		transformation.transform(transformedCoords[j]);
        	}
        	// get permutation of subunits and check validity/uniqueness
        	List<Integer> permutation = getPermutation();  
        	if (permutation.size() == 2 && isAllowedPermutation(permutation)) {
        		evaluateSolution(permutation, transformation, sphereAngle);
        	}
        }
    }
    
    private boolean isAllowedPermutation(List<Integer> permutation) {
    	List<Integer> seqClusterId = subunits.getSequenceClusterIds();
    	for (int i = 0; i < permutation.size(); i++) {
    		int j = permutation.get(i);
    		if (seqClusterId.get(i) != seqClusterId.get(j)) {
    			return false;
    		}
    	}
    	return true;
    }
    
    private void addEOperation() {
    	List<Integer> permutation = Arrays.asList(new Integer[]{0,1});
    	Matrix4d transformation = new Matrix4d();
    	transformation.setIdentity();
		combineWithTranslation(transformation);
    	AxisAngle4d axisAngle = new AxisAngle4d();
    	double rmsd = 0.0;
    	double caRmsd = 0.0;
    	int fold = 1; // ??
        Rotation rotation = createSymmetryOperation(permutation, transformation, axisAngle, rmsd, caRmsd, fold);
        rotations.addRotation(rotation);
    }
    
    /**
     * Transforms rotation axis so that it is perpendicular to the axis that connects the two subunit centers
     * 
     * @param axisAngle in x, y, z coordinate system
     * @return
     */
    private AxisAngle4d transformRotationAxis(AxisAngle4d axisAngle) {
    	Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
    	referenceTransformation.transform(v);
    	axisAngle.set(v, 0.0);
    	return axisAngle;
    }

    /**
     * Initializes a reference transformation that is used to transform the 
     * rotation axis into the subunit coordinate reference frame
     */
	private void initializeReferenceTransformation() {
		List<Point3d> centers = subunits.getCenters();
		
    	Vector3d subunitVector = new Vector3d(centers.get(0));
    	subunitVector.normalize();
    	
    	Vector3d orthogonalAxis = findOrthogonalAxis(subunitVector);
    	orthogonalAxis.normalize();

    	Point3d[] subunitReferenceFrame = {new Point3d(subunitVector), new Point3d(0.0, 0.0, 0.0), new Point3d(orthogonalAxis)};
    	Point3d[] cartesianReferenceFrame = {new Point3d(1.0, 0.0, 0.0), new Point3d(0.0, 0.0, 0.0), new Point3d(0.0, 0.0, 1.0)};
    	
    	AxisAngle4d a = new AxisAngle4d();
    	referenceTransformation = SuperPosition.superposeAtOrigin(cartesianReferenceFrame, subunitReferenceFrame, a);
	}

	/**
	 * Finds an orthogonal axis to the given vector
	 * @param v
	 * @return
	 */
	private Vector3d findOrthogonalAxis(Vector3d vector) {	
    	Vector3d orthogonalAxis = new Vector3d();
    	double dotProduct = 1.0;
    	
		MomentsOfInertia m = subunits.getMomentsOfInertia();
    	for (Vector3d axis: m.getPrincipalAxes()) {
    		if (Math.abs(vector.dot(axis)) < dotProduct) {
    			orthogonalAxis.cross(axis, vector);
    			dotProduct = Math.abs(vector.dot(axis));
    		}
    	}

		return orthogonalAxis;
	}

	private void evaluateSolution(List<Integer> permutation, Matrix4d transformation, AxisAngle4d axisAngle) {
		int fold = PermutationGroup.getOrder(permutation);
		double subunitRmsd = 0.0f; // should always be perfect for C2

		combineWithTranslation(transformation);
		double caRmsd = scorer.calcCalphaRMSD(transformation, permutation);
		if (caRmsd < 0.0 || caRmsd > parameters.getRmsdThreshold()) {
			return;
		}
		
		// if there is a better (lower RMSD) solution, remove the second solution and replace
		// it with this solution. Note, by convention, the first solution is E and should not be changed.
		if (caRmsd < bestRmsd && rotations.getOrder() == 2) {
			rotations.removeRotation(1);
			bestRmsd = caRmsd;
		}
		
		if (rotations.getOrder() == 1) {
			Rotation symmetryOperation = createSymmetryOperation(permutation, transformation, axisAngle, subunitRmsd, caRmsd, fold);
			rotations.addRotation(symmetryOperation);
		}
	}

    /**
     * Adds translational component to rotation matrix
     * @param rotTrans
     * @param rotation
     * @return
     */
    private void combineWithTranslation(Matrix4d rotation) {
        rotation.setTranslation(centroid);
        rotation.mul(rotation, centroidInverse);
    }

    private Rotation createSymmetryOperation(List<Integer> permutation, Matrix4d transformation, AxisAngle4d axisAngle, double subunitRmsd, double caRmsd, int fold) {
        Rotation s = new Rotation();
        s.setPermutation(new ArrayList<Integer>(permutation));
        s.setTransformation(new Matrix4d(transformation));
        s.setAxisAngle(new AxisAngle4d(axisAngle));
        s.setSubunitRmsd(subunitRmsd);
        s.setTraceRmsd(caRmsd);
        s.setFold(fold);
        return s;
    }

    private double calcDistanceThreshold() {
        double threshold = Double.MAX_VALUE;
        int n = subunits.getSubunitCount();
        List<Point3d> centers = subunits.getCenters();
        
        for (int i = 0; i < n - 1; i++) {
            Point3d pi = centers.get(i);
            for (int j = i + 1; j < n; j++) {
                Point3d pj = centers.get(j);
                threshold = Math.min(threshold, pi.distanceSquared(pj));
            }
        }
        double distanceThreshold = Math.sqrt(threshold);

//        System.out.println("Distance threshold: " + distanceThreshold);
        distanceThreshold = Math.max(distanceThreshold, parameters.getRmsdThreshold());
        
        return distanceThreshold;
    }

    private List<Integer> getPermutation() {
        List<Integer> permutation = new ArrayList<Integer>(transformedCoords.length);
        double sum = 0.0f;

        for (Point3d t: transformedCoords) {
            int closest = -1;
            double minDist = Double.MAX_VALUE;

            for (int j = 0; j< originalCoords.length; j++) {
            	double dist = t.distanceSquared(originalCoords[j]);
                if (dist < minDist) {
                    closest = j;
                    minDist = dist;
                }
            }
            
            sum += minDist;
            if (closest == -1) {
         	   break;
            }
            permutation.add(closest);
        }
        double rmsd = Math.sqrt(sum / transformedCoords.length);

        if (rmsd > distanceThreshold || permutation.size() != transformedCoords.length) {
            permutation.clear();
            return permutation;
        }

        // check uniqueness of indices
        Set<Integer> set = new HashSet<Integer>(permutation);
        
        // if size mismatch, clear permutation (its invalid)
        if (set.size() != originalCoords.length) {
            permutation.clear();
        }
      
        return permutation;
    }

    private void initialize() {
        scorer = new QuatSuperpositionScorer(subunits);
        
        // translation to centered coordinate system
        centroid = new Vector3d(subunits.getCentroid());

        // translation back to original coordinate system
        Vector3d reverse = new Vector3d(centroid);
        reverse.negate();
        centroidInverse.set(reverse);
        // On LINUX there seems to be a bug with vecmath, and element m33 is zero. Here we make sure it's 1.
        centroidInverse.setElement(3, 3, 1.0);

        List<Point3d> centers = subunits.getCenters();
        int n = subunits.getSubunitCount();

        originalCoords = new Point3d[n];
        transformedCoords = new Point3d[n];

        for (int i = 0; i < n; i++) {
            originalCoords[i] = centers.get(i);
            transformedCoords[i] = new Point3d();
        }

        distanceThreshold = calcDistanceThreshold();
        initializeReferenceTransformation();
    }
}
