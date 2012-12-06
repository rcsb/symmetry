
package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.quaternary.geometry.DistanceBox;
import org.biojava3.structure.quaternary.geometry.MomentsOfInertia;
import org.biojava3.structure.quaternary.geometry.SphereSampler;
import org.biojava3.structure.quaternary.geometry.SuperPosition;


/**
 *
 * @author Peter
 */
public class RotationSolver implements QuatSymmetrySolver {
    private Subunits subunits = null;
    private QuatSymmetryParameters parameters = null;
    
 //   private double rmsdThreshold = 0.0f;
    private double distanceThreshold = 0.0f;
    private DistanceBox<Integer> box = null;
    private Vector3d centroid = new Vector3d();
    private Matrix4d centroidInverse = new Matrix4d();
    private Point3d[] originalCoords = null;
    private Point3d[] transformedCoords = null;
    private Set<List<Integer>> hashCodes = new HashSet<List<Integer>>();

    private RotationGroup rotations = new RotationGroup();
    private QuatSuperpositionScorer scorer = null;
    
//    private SuperPositionQCP sp = null;

    public RotationSolver(Subunits subunits, QuatSymmetryParameters parameters) {
    	if (subunits.getSubunitCount()== 2) {
    		throw new IllegalArgumentException("RotationSolver cannot be applied to subunits with 2 centers");
    	}
        this.subunits = subunits;
        this.parameters = parameters;
//        sp = new SuperPositionQCP();
    }

	public RotationGroup getSymmetryOperations() {
        if (rotations.getOrder() == 0) {
            solve();
            completeRotationGroup();
        }
        return rotations;
    }

    private void solve() {
        initialize();
        
        int maxSymOps = subunits.getSubunitCount();
        // for cases with icosahedral symmetry n cannot be higher than 60, should check for spherical symmetry here
        // isSpherical check added 08-04-11
        if (maxSymOps % 60 == 0 && isSpherical()) {
            maxSymOps = 60;
         }

        AxisAngle4d sphereAngle = new AxisAngle4d();
        Matrix4d transformation = new Matrix4d();

        int n = subunits.getSubunitCount();
        List<Double> angles = getAngles();

       for (int i = 0; i < SphereSampler.getSphereCount(); i++) {
            SphereSampler.getAxisAngle(i, sphereAngle);
            for (double angle : angles) {
                // apply rotation
                sphereAngle.angle = angle;
                transformation.set(sphereAngle);
                for (int j = 0; j < n; j++) {
                    transformedCoords[j].set(originalCoords[j]);
                    transformation.transform(transformedCoords[j]);
                }

                // get permutation of subunits and check validity/uniqueness             
                List<Integer> permutation = getPermutation();
  //              System.out.println("Rotation Solver: permutation: " + i + ": " + permutation);
                if (! isValidPermutation(permutation)) {
                    continue;
                }
               
                boolean newPermutation = evaluatePermutation(permutation);
                if (newPermutation) {
                	completeRotationGroup();
                }
                
                // check if all symmetry operations have been found.          
                if (rotations.getOrder() >= maxSymOps) {
                	return;
                }
            }
        }
    }
    
    private void completeRotationGroup() {
    	PermutationGroup g = new PermutationGroup();
    	for (int i = 0; i < rotations.getOrder(); i++) {
    		Rotation s = rotations.getRotation(i);
    		g.addPermutation(s.getPermutation());
    	}
    	g.completeGroup();
    	
    	// the group is complete, nothing to do
    	if (g.getOrder() == rotations.getOrder()) {
    		return;
    	}
    	
 //   	System.out.println("complete group: " +  rotations.getOrder() +"/" + g.getOrder());
    	// try to complete the group
    	for (int i = 0; i < g.getOrder(); i++) {
    		List<Integer> permutation = g.getPermutation(i);
    		if (isValidPermutation(permutation)) {
    			  // perform permutation of subunits
                evaluatePermutation(permutation);
    		}
    	}
    }

	private boolean evaluatePermutation(List<Integer> permutation) {
		// permutate subunits
		for (int j = 0, n = subunits.getSubunitCount(); j < n; j++) {
		    transformedCoords[j].set(originalCoords[permutation.get(j)]);
		}

		int fold = PermutationGroup.getOrder(permutation);
		// get optimal transformation and axisangle by superimposing subunits
		AxisAngle4d axisAngle = new AxisAngle4d();
		
		// ---
//		long t1 = System.nanoTime();
//		// are these coordinates precentered??
//		sp.set(transformedCoords, originalCoords);
//		sp.setCentered(true);
//		double subunitRmsd = sp.getRmsd();
//		long t2 = System.nanoTime();
	    // ----
		Matrix4d transformation = SuperPosition.superposeAtOrigin(transformedCoords, originalCoords, axisAngle);
		double subunitRmsd = SuperPosition.rmsd(transformedCoords, originalCoords);
		//
	
		
	//	long t3 = System.nanoTime();
		// --
//		System.out.println(" new rmsd: " + subunitRmsd);
//		System.out.println(" new: " + (t2-t1));
		// --
		
		if (subunitRmsd < parameters.getRmsdThreshold()) {
//			Matrix4d transformation = sp.getTransformationMatrix();
//			transformedCoords = sp.getTransformedCoordinates();
//			axisAngle.set(transformation);
//			
//			long t3 = System.nanoTime();
//			System.out.println(" total time: " + (t3-t1));
			// transform to original coordinate system
			combineWithTranslation(transformation);
			// evaluate superposition of CA traces with GTS score
			double caRmsd = scorer.calcCalphaRMSD(transformation, permutation);
			if (caRmsd < 0.0 || caRmsd > parameters.getRmsdThreshold()) {
				return false;
			}
			Rotation symmetryOperation = createSymmetryOperation(permutation, transformation, axisAngle, subunitRmsd, caRmsd, fold);
			rotations.addRotation(symmetryOperation);
			return true;
		}
		return false;
	}

    private List<Double> getAngles() {
        int n = subunits.getSubunitCount();
        // for spherical symmetric cases, n cannot be higher than 60
        if (n % 60 == 0 && isSpherical()) {
           n = 60;
        }
        List<Integer> folds = subunits.getFolds();
        List<Double> angles = new ArrayList<Double>(folds.size()-1);

        // note this loop starts at 1, we do ignore 1-fold symmetry, which is the first entry
        for (int fold: folds) {
        	if (fold > 0 && fold <= n) {
        		angles.add(2* Math.PI/fold);
        	}
        }
        return angles;
    }
    
    private boolean isSpherical() {
    	MomentsOfInertia m = subunits.getMomentsOfInertia();
    	return m.getSymmetryClass(0.05) == MomentsOfInertia.SymmetryClass.SYMMETRIC;
    }

    private boolean isValidPermutation(List<Integer> permutation) {
    	  // if this permutation is a duplicate, return false
    	if (hashCodes.contains(permutation)) {
    		return false;
    	}
        if (permutation.size() == 0) {
 //       	System.out.println("permutation size zero");
            return false;
        }
        // check if permutation is allowed
        if (! isAllowedPermutation(permutation)) {
        	return false;
        }
     // get fold and make sure there is only one E (fold=1) permutation
        int fold = PermutationGroup.getOrder(permutation);
        if (rotations.getOrder() > 1 && fold == 1) {
//        	System.out.println("Symop = 1");
            return false;
        }
        if (fold == 0 || subunits.getSubunitCount() % fold != 0) {
  //      	System.out.println(permutation);
  //      	System.out.println("Remove: " + subunits.getSubunitCount() + " / " + fold);
        	return false;
        }
        
        // if this permutation is a duplicate, returns false
        return hashCodes.add(permutation);
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

    private Rotation createSymmetryOperation(List<Integer> permutation, Matrix4d transformation, AxisAngle4d axisAngle, double subunitRmsd, double rmsd, int fold) {
        Rotation s = new Rotation();
        s.setPermutation(new ArrayList<Integer>(permutation));
        s.setTransformation(new Matrix4d(transformation));
        s.setAxisAngle(new AxisAngle4d(axisAngle));
        s.setSubunitRmsd(subunitRmsd);
        s.setTraceRmsd(rmsd);
        s.setFold(fold);
        return s;
    }

    private void setupDistanceBox() {
        distanceThreshold = calcDistanceThreshold();
        box = new DistanceBox<Integer>(distanceThreshold);

        for (int i = 0; i < originalCoords.length; i++) {
            box.addPoint(originalCoords[i], i);
        }
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

 //       System.out.println("Distance threshold: " + distanceThreshold);
        distanceThreshold = Math.max(distanceThreshold, parameters.getRmsdThreshold());
        
        return distanceThreshold;
    }

    private List<Integer> getPermutation() {
        List<Integer> permutation = new ArrayList<Integer>(transformedCoords.length);
        double sum = 0.0f;

        for (Point3d t: transformedCoords) {
            List<Integer> neighbors = box.getNeighborsWithCache(t);
            int closest = -1;
            double minDist = Double.MAX_VALUE;

           for (int j : neighbors) {
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
  //      	System.out.println("RotationSolver: getPermutation: duplicate members" + set.size());
            permutation.clear();
        }
      
//        System.out.println("RMSD: " + rmsd + " missed: " + Math.sqrt(missed) + " closest: " + missedI + " count: " + count);
 //       System.out.println("P1: " + permutation);
 //       System.out.println("P2: " + p2);
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

        List<Point3d> centers = subunits.getCenters();
        int n = subunits.getSubunitCount();

        originalCoords = new Point3d[n];
        transformedCoords = new Point3d[n];

        for (int i = 0; i < n; i++) {
            originalCoords[i] = centers.get(i);
            transformedCoords[i] = new Point3d();
        }

        setupDistanceBox();
    }
}
