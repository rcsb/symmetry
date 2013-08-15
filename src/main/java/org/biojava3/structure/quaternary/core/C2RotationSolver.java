
package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.quaternary.geometry.SuperPosition;

/**
 *
 * @author Peter
 */
public class C2RotationSolver implements QuatSymmetrySolver {
    private Subunits subunits = null;
    private QuatSymmetryParameters parameters = null;
    private Vector3d centroid = new Vector3d();
    private Matrix4d centroidInverse = new Matrix4d();
    private QuatSuperpositionScorer scorer = null;

    private RotationGroup rotations = new RotationGroup();


    public C2RotationSolver(Subunits subunits, QuatSymmetryParameters parameters) {
        if (subunits.getSubunitCount() != 2) {
    		throw new IllegalArgumentException("C2RotationSolver can only be applied to cases with 2 centers");
    	}
        this.subunits = subunits;
        this.parameters = parameters;
        this.scorer = new QuatSuperpositionScorer(subunits);
    }
   
    public RotationGroup getSymmetryOperations() {
        if (rotations.getOrder() == 0) {
            solve();
        }
        return rotations;
    }
    
    private void solve() {
    	initialize();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		List<Point3d[]> traces = subunits.getTraces();

//		Point3d[] x = SuperPosition.clonePoint3dArray(traces.get(0));
//		SuperPosition.center(x);
//		Point3d[] y = SuperPosition.clonePoint3dArray(traces.get(1));
//		SuperPosition.center(y);
		
		Point3d[] x = SuperPosition.clonePoint3dArray(traces.get(0));
		SuperPosition.translate(new Point3d(trans), x);
		Point3d[] y = SuperPosition.clonePoint3dArray(traces.get(1));
		SuperPosition.translate(new Point3d(trans), y);

		AxisAngle4d axisAngle = new AxisAngle4d();

		Matrix4d transformation = SuperPosition.superposeAtOrigin(x, y, axisAngle);
		double caRmsd = SuperPosition.rmsd(x,  y);
		
		// TODO this is not the proper way to calculate the TM score. Each subunit should be treated
		// separately!
		double caTmScoreMin = SuperPosition.TMScore(x, y, x.length);
		
		// if rmsd or angle deviation is above threshold, stop
		double angleThresholdRadians = Math.toRadians(parameters.getAngleThreshold());
		double deltaAngle = Math.abs(Math.PI-axisAngle.angle);
	
		if (caRmsd > parameters.getRmsdThreshold() || deltaAngle > angleThresholdRadians) {
			rotations.setC1(subunits.getSubunitCount());
			return;
		}
		
		// add unit operation
		addEOperation();

		// add C2 operation
		List<Integer> permutation = new ArrayList<Integer>();
		permutation.add(new Integer(1));
		permutation.add(new Integer(0));

		// combine with translation
		Matrix4d combined = new Matrix4d();
		combined.setIdentity();
		combined.setTranslation(trans);
		transformation.mul(combined);
	
		Rotation symmetryOperation = createSymmetryOperation(permutation, transformation, axisAngle, 0.0, caRmsd, caTmScoreMin, 2);
		rotations.addRotation(symmetryOperation);
    }
    
    private void addEOperation() {
    	List<Integer> permutation = Arrays.asList(new Integer[]{0,1});
    	Matrix4d transformation = new Matrix4d();
    	transformation.setIdentity();
		combineWithTranslation(transformation);
    	AxisAngle4d axisAngle = new AxisAngle4d();
    	double rmsd = 0.0;
    	double caRmsd = 0.0;
    	double caTmScoreMin = 1.0;
    	int fold = 1; // ??
        Rotation rotation = createSymmetryOperation(permutation, transformation, axisAngle, rmsd, caRmsd, caTmScoreMin, fold);
        rotations.addRotation(rotation);
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

    private Rotation createSymmetryOperation(List<Integer> permutation, Matrix4d transformation, AxisAngle4d axisAngle, double subunitRmsd, double caRmsd, double caTmScoreMin, int fold) {
        Rotation s = new Rotation();
        s.setPermutation(new ArrayList<Integer>(permutation));
        s.setTransformation(new Matrix4d(transformation));
        s.setAxisAngle(new AxisAngle4d(axisAngle));
        s.setSubunitRmsd(subunitRmsd);
        s.setTraceRmsd(caRmsd);
        s.setTraceTmScoreMin(caTmScoreMin);
        s.setFold(fold);
        return s;
    }

    private void initialize() {     
        // translation to centered coordinate system
        centroid = new Vector3d(subunits.getCentroid());
       // translation back to original coordinate system
        Vector3d reverse = new Vector3d(centroid);
        reverse.negate();     
        centroidInverse.set(reverse);
//        // On LINUX there seems to be a bug with vecmath, and element m33 is zero. Here we make sure it's 1.
        centroidInverse.setElement(3, 3, 1.0);
    }

}
