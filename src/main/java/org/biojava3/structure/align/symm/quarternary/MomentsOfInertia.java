
package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.List;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.jama.EigenvalueDecomposition;
import org.biojava.bio.structure.jama.Matrix;

/**
 *
 * @author Peter
 */
public class MomentsOfInertia {
    private List<Point3d> points = new ArrayList<Point3d>();
    private List<Double> masses = new ArrayList<Double>();
    
    private boolean modified = true;
    
    private double[] principalMomentsOfInertia = new double[3];
    private Vector3d[] principalAxes = new Vector3d[3];

    public enum SymmetryClass {LINEAR, PROLATE, OBLATE, SYMMETRIC, ASYMMETRIC};
    
    /** Creates a new instance of MomentsOfInertia */
    public MomentsOfInertia() {
    }
    
    public void addPoint(Point3d point, double mass) {
        points.add(point);
        masses.add(mass);
        modified = true;
    }
    
    public Point3d centerOfMass() {
        Point3d center = new Point3d();
        
        double totalMass = 0.0f;
        for (int i = 0, n = points.size(); i < n; i++) {
            double mass = masses.get(i);
            totalMass += mass;
            center.scaleAdd(mass, points.get(i), center);
        }
        center.scale(1.0f/totalMass);
        return center;
    }

    public Point3d centerOfMassD() {
        Point3d center = new Point3d();
        Point3d temp = new Point3d();

        double totalMass = 0.0f;
        for (int i = 0, n = points.size(); i < n; i++) {
            double mass = masses.get(i);
            totalMass += mass;
            temp.set(points.get(i));
            center.scaleAdd(mass, temp, center);
        }
        center.scale(1.0f/totalMass);
        return center;
    }
    
    public double[] getPrincipalMomentsOfInertia() {
        if (modified) {
            diagonalizeTensor();
            modified = false;
        }
        return principalMomentsOfInertia;
    }
    
    public Vector3d[] getPrincipalAxes() {
        if (modified) {
            diagonalizeTensor();
            modified = false;
        }
        return principalAxes;
    }

    public SymmetryClass getSymmetryClass(double threshold) {
        if (modified) {
            diagonalizeTensor();
            modified = false;
        }
        double ia = principalMomentsOfInertia[0];
        double ib = principalMomentsOfInertia[1];
        double ic = principalMomentsOfInertia[2];
        boolean c1 = (ib - ia) / (ib + ia) < threshold;
        boolean c2 = (ic - ib) / (ic + ib) < threshold;

        if (c1 && c2) {
            return SymmetryClass.SYMMETRIC;
        }
        if (c1) {
        	return SymmetryClass.OBLATE;
        }
        if (c2) {
        	return SymmetryClass.PROLATE;
        }
        return SymmetryClass.ASYMMETRIC;
    }

    public double symmetryCoefficient() {
        if (modified) {
            diagonalizeTensor();
            modified = false;
        }
        double ia = principalMomentsOfInertia[0];
        double ib = principalMomentsOfInertia[1];
        double ic = principalMomentsOfInertia[2];
        double c1 = 1.0f - (ib - ia) / (ib + ia);
        double c2 = 1.0f - (ic - ib) / (ic + ib);
        return Math.max(c1, c2);
//        return Math.min(c1, c2);
//        return 0.5 * (c1 + c2);
    }
    
    public double getAsymmetryParameter(double threshold) {
       if (modified) {
            diagonalizeTensor();
            modified = false;
        }
       if (getSymmetryClass(threshold).equals(SymmetryClass.SYMMETRIC)) {
    	   return 0.0;
//           throw new IllegalStateException("Asymmetry parameter is undefined for a symmetric top.");
       }
       double a = 1.0/principalMomentsOfInertia[0];
       double b = 1.0/principalMomentsOfInertia[1];
       double c = 1.0/principalMomentsOfInertia[2];
       return (2 * b - a - c) / (a - c);
    }
    
    private double[][] getInertiaTensor() {
//        http://en.wikipedia.org/wiki/Moment_of_inertia
        Point3d p = new Point3d();
        
        double[][] tensor = new double[3][3];
        
        // calculate the inertia tensor at center of mass
        Point3d com = centerOfMassD();
        
        for (int i = 0, n = points.size(); i < n; i++) {
            double mass = masses.get(i);
            p.sub(points.get(i), com);
            double px = p.x;
            double py = p.y;
            double pz = p.z;
            
            tensor[0][0] += mass * (py * py + pz * pz);
            tensor[1][1] += mass * (px * px + pz * pz);
            tensor[2][2] += mass * (px * px + py * py);
            
            tensor[0][1] -= mass * px * py;
            tensor[0][2] -= mass * px * pz;
            tensor[1][2] -= mass * py * pz;
        }
        
        tensor[1][0] = tensor[0][1];
        tensor[2][0] = tensor[0][2];
        tensor[2][1] = tensor[1][2];
        
        return tensor;
    }
    
    private void diagonalizeTensor() {
        Matrix m = new Matrix(getInertiaTensor());
        
        EigenvalueDecomposition eig = m.eig();
        double[] eigenValues = eig.getRealEigenvalues();
        double[][] eigenVectors = eig.getV().getArray();
        
        for (int i = 0; i < 3; i++) {
            principalMomentsOfInertia[i] = eigenValues[i];
            
            principalAxes[i] = new Vector3d(
                    eigenVectors[i][0],
                    eigenVectors[i][1],
                    eigenVectors[i][2]);
        }
    }
}
