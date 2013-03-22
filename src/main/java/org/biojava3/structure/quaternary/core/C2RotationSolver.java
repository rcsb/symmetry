
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
public class C2RotationSolver implements QuatSymmetrySolver {
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

    public C2RotationSolver(Subunits subunits, QuatSymmetryParameters parameters) {
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
        initialize();
        
        // add the unit operation
        addEOperation();
        
        AxisAngle4d sphereAngle = new AxisAngle4d();
        
        Matrix4d transformation = new Matrix4d();
        transformation.setIdentity();
        
        int n = subunits.getSubunitCount();
        System.out.println("C2 rotation solver subunits: " + n + " " + subunits.getSequenceClusterIds());
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
        
        System.out.println("C2RotationSolver #rotations: " + rotations.getOrder());
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
    	Point3d[] cartesianReferenceFrame = {new Point3d(1.0, 0.0, 0.0), new Point3d(0.0, 1.0, 0.0), new Point3d(0.0, 0.0, 1.0)};
    	
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
		
		//System.out.println("C2RotationSolver evaluateSolution: " + fold + " " + permutation + " " + transformation.getElement(3,3));
		
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
		
		System.out.println(String.format("C2RotationSolver evaluateSolution: %10.2f %10.2f %d", caRmsd , bestRmsd , rotations.getOrder()));
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
