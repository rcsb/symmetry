
package org.biojava3.structure.align.symm.quarternary;

import java.util.Arrays;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.jama.EigenvalueDecomposition;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.jama.SingularValueDecomposition;

/**
 *
 * @author Peter
 */
public final class SuperPosition {
    
    // this class cannot be instantiated
    private SuperPosition() {}
    
    
    public static Matrix4d superposeSVD(Point3d[] x, Point3d[] y) {
    	Atom[] a = new Atom[x.length];
    	Atom[] b = new Atom[x.length];
    	System.out.println("Superpose:");
    	System.out.println(Arrays.toString(x));
    	System.out.println(Arrays.toString(y));
    	
    	for (int i = 0; i < x.length; i++) {
    		double[] c = new double[3];
    		x[i].get(c);
    		Atom ax = new AtomImpl();
    		ax.setCoords(c);
    		a[i] = ax;
    		y[i].get(c);
    		Atom bx = new AtomImpl();
    		bx.setCoords(c);
    		b[i] = bx;
    	}
    	SVDSuperimposer s = null;
    	try {
		     s = new SVDSuperimposer(a, b);
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Matrix rot = s.getRotation();
		Atom t = s.getTranslation();
		Matrix4d rotTrans = new Matrix4d();
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				rotTrans.setElement(i, j, rot.get(i, j));
			}
		}
//		System.out.println(rotTrans);
		rotTrans.setTranslation(new Vector3d(t.getX(), t.getY(), t.getZ()));	
		 // transform coordinates
        transform(rotTrans, x);
		System.out.println(rotTrans);

        return rotTrans;
    }
    
    public static Matrix4d superpose(Point3d[] x, Point3d[] y) {
        //superpose x onto y
        Point3d[] ref = clonePoint3dArray(y);
        
        Point3d ytrans = centroid(ref);
        ytrans.negate();
        translate(ytrans, ref);
        
        center(x);
        
        // calculate quaternion from relative orientation
        Quat4d q = quaternionOrientation(x, ref);
        q.conjugate();
        
        Matrix4d rotTrans = new Matrix4d();
        rotTrans.set(q);
        
        // set translational component of transformation matrix
        ytrans.negate();
        rotTrans.setTranslation(new Vector3d(ytrans));
        
        // tranform coordinates
        transform(rotTrans, x);

        return rotTrans;
    }

    public static Matrix4d superposeAtOrigin(Point3d[] x, Point3d[] y, AxisAngle4d axisAngle) {
        Quat4d q = quaternionOrientation(x, y);
        q.conjugate();
        
        Matrix4d rotTrans = new Matrix4d();
        rotTrans.set(q);
        axisAngle.set(q);
        transform(rotTrans, x);
//		System.out.println(rotTrans);
        
        return rotTrans;
    }
    
    public static double rmsd(Point3d[] x, Point3d[] y) {
        double sum = 0.0f;
        for (int i = 0; i < x.length; i++) {
            sum += x[i].distanceSquared(y[i]);
        }
        return (double)Math.sqrt(sum/x.length);
    }

    public static double rmsdMin(Point3d[] x, Point3d[] y) {
        double sum = 0.0f;
        for (int i = 0; i < x.length; i++) {
            double minDist = Double.MAX_VALUE;
            for (int j = 0; j < y.length; j++) {
               minDist = Math.min(minDist, x[i].distanceSquared(y[j]));
            }
            sum += minDist;
        }
        return (double)Math.sqrt(sum/x.length);
    }

    public static double GTSlikeScore(Point3d[] x, Point3d[] y) {
        int contacts = 0;

        for (Point3d px: x) {
            double minDist = Double.MAX_VALUE;
            
            for (Point3d py: y) {
               minDist = Math.min(minDist, px.distanceSquared(py));
            }
            
            if (minDist > 64) continue;
            contacts++;

            if (minDist > 16) continue;
            contacts++;

            if (minDist > 4) continue;
            contacts++;

            if (minDist > 1) continue;
            contacts++;
        }

        return contacts*25/x.length;
    }
    
    public static int contacts(Point3d[] x, Point3d[] y, double maxDistance) {
        int contacts = 0;
        for (int i = 0; i < x.length; i++) {
            double minDist = Double.MAX_VALUE;
            for (int j = 0; j < y.length; j++) {
               minDist = Math.min(minDist, x[i].distanceSquared(y[j]));
            }
            if (minDist < maxDistance*maxDistance) {
                contacts++;
            }
        }
        return contacts;
    }
    
    public static void transform(Matrix4d rotTrans, Point3d[] x) {
        for (Point3d p: x) {
            rotTrans.transform(p);
        }
    }
    
    public static void translate(Point3d trans, Point3d[] x) {
        for (Point3d p: x) {
            p.add(trans);
        }
    }
    
    public static void center(Point3d[] x) {
        Point3d center = centroid(x);
        center.negate();
        translate(center, x);
    }
    
    public static Point3d centroid(Point3d[] x) {
        Point3d center = new Point3d();
        for (Point3d p: x) {
            center.add(p);
        }
        center.scale(1.0f/x.length);
        return center;
    }
    
    private static Quat4d quaternionOrientationSVD(Point3d[] a, Point3d[] b)  {
        Matrix m = calcFormMatrix(a, b);
        SingularValueDecomposition eig = m.svd();
        double[][] v = eig.getV().getArray();
//        System.out.println(eig.getV());   
//        System.out.println("q: " + v[0][3] + " " + v[1][3] +" " + v[2][3] + " " + v[3][3]);
        Quat4d q = new Quat4d(v[2][3], v[1][3], v[0][3], v[3][3]);
 //       System.out.println("q: " + q); 
        return q;
    }
    
    private static Quat4d quaternionOrientation(Point3d[] a, Point3d[] b)  {
        Matrix m = calcFormMatrix(a, b);
        EigenvalueDecomposition eig = m.eig();
//        double[] e = eig.getRealEigenvalues();
//        System.out.println("eigval: " + Arrays.toString(e));
        double[][] v = eig.getV().getArray();
//        System.out.println(eig.getV()); 
        Quat4d q = new Quat4d(v[1][3], v[2][3], v[3][3], v[0][3]);
//        System.out.println("q: " + q); 
        return q;
    }
    
    private static Matrix calcFormMatrix(Point3d[] a, Point3d[] b) {
        double xx=0.0, xy=0.0, xz=0.0, yx=0.0, yy=0.0, yz=0.0, zx=0.0, zy=0.0, zz=0.0;
        
        for (int i = 0; i < a.length; i++) {
            xx += a[i].x * b[i].x;
            xy += a[i].x * b[i].y;
            xz += a[i].x * b[i].z;
            yx += a[i].y * b[i].x;
            yy += a[i].y * b[i].y;
            yz += a[i].y * b[i].z;
            zx += a[i].z * b[i].x;
            zy += a[i].z * b[i].y;
            zz += a[i].z * b[i].z;
        }
        
        double[][] f = new double[4][4];
        f[0][0] = xx + yy + zz;
        f[0][1] = zy - yz;
        f[1][0] = f[0][1];
        f[1][1] = xx - yy - zz;
        f[0][2] = xz - zx;
        f[2][0] = f[0][2];
        f[1][2] = xy + yx;
        f[2][1] = f[1][2];
        f[2][2] = yy - zz - xx;
        f[0][3] = yx - xy;
        f[3][0] = f[0][3];
        f[1][3] = zx + xz;
        f[3][1] = f[1][3];
        f[2][3] = yz + zy;
        f[3][2] = f[2][3];
        f[3][3] = zz - xx - yy;
        
        return new Matrix(f);
    }
    
    public static Point3d[] clonePoint3dArray(Point3d[] x) {
        Point3d[] clone = new Point3d[x.length];
        for (int i = 0; i < x.length; i++) {
           clone[i] = new Point3d(x[i]);
        }
        return clone;
    }
  
}
