
package org.biojava3.structure.align.symm.quaternary;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;


/**
 *
 * @author Peter
 */
public class SystematicSolver implements QuatSymmetrySolver {
    private Subunits subunits = null;
//    private double subunitRmsdThreshold = 3.0f;
    private double rmsdThreshold = 6.0f;
    private Point3d[] originalCoords = null;
    private Point3d[] transformedCoords = null;
    private RotationGroup rotations = new RotationGroup();
    private Vector3d centroid = new Vector3d();
    private Matrix4d centroidInverse = new Matrix4d();
    private QuatSuperpositionScorer scorer = null;
    private Set<List<Integer>> hashCodes = new HashSet<List<Integer>>();

    public SystematicSolver(Subunits subunits, double rmsdThreshold) {
    	if (subunits.getSubunitCount()== 2) {
    		throw new IllegalArgumentException("SystematicSolver cannot be applied to subunits with 2 centers");
    	}
        this.subunits = subunits;
        this.rmsdThreshold = rmsdThreshold;
    }

//    public void setRmsdThreshold(double rmsdThreshold) {
//        this.rmsdThreshold = rmsdThreshold;
//    }

    public RotationGroup getSymmetryOperations() {
        if (rotations.getOrder() == 0) {
            solve();
        }
        return rotations;
    }

    private void solve() {
        initialize();
        int n = subunits.getSubunitCount();
        PermutationGenerator g = new PermutationGenerator(n);

        // loop over all permutations
        while (g.hasMore()) {
            int[] perm = g.getNext();
            List<Integer> permutation = new ArrayList<Integer>(perm.length);
            for (int j = 0; j < n; j++) {
                permutation.add(perm[j]);
            }
            
            if (! isValidPermutation(permutation)) {
                continue;
            }
            
            boolean newPermutation = evaluatePermutation(permutation);
            if (newPermutation) {
            	completeRotationGroup();
            }
            
            if (rotations.getOrder() >= subunits.getSubunitCount()) {
            	return;
            }
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

    private Rotation createSymmetryOperation(List<Integer> permutation, Matrix4d transformation, AxisAngle4d axisAngle, double rmsd, double gts, double caRmsd, int fold) {
        Rotation s = new Rotation();
        s.setPermutation(new ArrayList<Integer>(permutation));
        s.setTransformation(new Matrix4d(transformation));
        s.setAxisAngle(new AxisAngle4d(axisAngle));
        s.setSubunitRmsd(rmsd);
        s.setTraceRmsd(caRmsd);
        s.setFold(fold);
        return s;
    }
    
    private void completeRotationGroup() {
    	PermutationGroup g = new PermutationGroup();
    	for (int i = 0; i < rotations.getOrder(); i++) {
    		Rotation s = rotations.getRotation(i);
    		g.addPermutation(s.getPermutation());
    	}
    	g.completeGroup();
    	
 //   	System.out.println("Completing rotation group from: " +symmetryOperations.getSymmetryOperationCount() + " to " + g.getPermutationCount());
    	
    	// the group is complete, nothing to do
    	if (g.getOrder() == rotations.getOrder()) {
    		return;
    	}
    	
  //  	System.out.println("complete group: " +  rotations.getOrder() +"/" + g.getOrder());
    	// try to complete the group
    	for (int i = 0; i < g.getOrder(); i++) {
    		List<Integer> permutation = g.getPermutation(i);
    		if (isValidPermutation(permutation)) {
    			  // perform permutation of subunits
                evaluatePermutation(permutation);
    		}
    	}
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
        
        // check if permutation is pseudosymmetric
        if (! isAllowedPermuation(permutation)) {
        	return false;
        }
        
        // get fold and make sure there is only one E (fold=1) permutation
        int fold = PermutationGroup.getOrder(permutation);
        if (rotations.getOrder() > 1 && fold == 1) {
            return false;
        }
        if (fold == 0 || subunits.getSubunitCount() % fold != 0) {
 //       	System.out.println("Remove: " + subunits.getSubunitCount() + " / " + fold);
        	return false;
        }
        
        // if this permutation is a duplicate, returns false
        return hashCodes.add(permutation);
    }

    private boolean isAllowedPermuation(List<Integer> permutation) {
    	List<Integer> seqClusterId = subunits.getSequenceClusterIds();
    	for (int i = 0; i < permutation.size(); i++) {
    		int j = permutation.get(i);
    		if (seqClusterId.get(i) != seqClusterId.get(j)) {
    			return false;
    		}
    	}
    	return true;
    }
    
	private boolean evaluatePermutation(List<Integer> permutation) {
		// permutate subunits
		for (int j = 0, n = subunits.getSubunitCount(); j < n; j++) {
		    transformedCoords[j].set(originalCoords[permutation.get(j)]);
		}

		int fold = PermutationGroup.getOrder(permutation);
		// get optimal transformation and axisangle by superimposing subunits
		AxisAngle4d axisAngle = new AxisAngle4d();
		Matrix4d transformation = SuperPosition.superposeAtOrigin(transformedCoords, originalCoords, axisAngle);
		double rmsd = SuperPosition.rmsd(transformedCoords, originalCoords);
		// handle ambiguity of superposing two subunits with the E operation (count=1)
        if (subunits.getSubunitCount() == 2 && rotations.getOrder() == 0) {
        	transformation.setIdentity();
        }
 //               System.out.println("Complete: " + permutation + " rmsd: " + rmsd);
		// check if it meets criteria and save symmetry operation
		if (rmsd < rmsdThreshold) {
			// transform to original coordinate system
		    combineWithTranslation(transformation);
		    // evaluate superposition of CA traces with GTS score
//		    double gts = scorer.calcGtsMinScore(transformation, permutation);
		    double gts = 0;
//                    System.out.println("Complete: " + permutation + " gts: " + gts);
//		    if (gts > gtsThreshold) {
		    	double caRmsd = scorer.calcCalphaRMSD(transformation, permutation);
		    	if (caRmsd < 0.0) {
		    		return false;
		    	}
		    	if (caRmsd > rmsdThreshold) {
		            return false;
		    	}
		        Rotation symmetryOperation = createSymmetryOperation(permutation, transformation, axisAngle, rmsd, gts, caRmsd, fold);
		        rotations.addRotation(symmetryOperation);
		        return true;
//		    }
		}
		return false;
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
    }
}
