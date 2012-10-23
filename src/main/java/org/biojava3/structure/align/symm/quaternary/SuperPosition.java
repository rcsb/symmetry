
package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.jama.EigenvalueDecomposition;
import org.biojava.bio.structure.jama.Matrix;

/**
 *
 * @author Peter
 */
public final class SuperPosition {
    
    private Point3d[] a = null;
    private Point3d[] b = null;
    private double[] weight = null;
    private boolean rmsdOnly = true;
    private double maxRmsd = Double.MAX_VALUE;
    private Matrix3d rotmat = new Matrix3d();
    private Matrix4d transformation = new Matrix4d();
    private double rmsd = 0;
    
   
    public SuperPosition(boolean rmsdOnly, double maxRmsd) {
    	this.rmsdOnly = rmsdOnly;
    	this.maxRmsd = maxRmsd;
    }
    
    public void superposeNew(Point3d[] a, Point3d[] b) {
    	this.a = a;
    	this.b = b;
    	calcRMSDRotationalMatrix(a, b);
    }
    
    public void superposeNew(Point3d[] a, Point3d[] b, double[] mass) {
    	this.a = a;
    	this.b = b;
    	this.weight = mass;
    	calcRMSDRotationalMatrix(a, b);
    }
    
    public double getRmsd() {
    	return rmsd;
    }
    
    public Matrix4d getTransformationMatrix() {
    	return transformation; 	
    }
    public Matrix3d getRotationMatrix() {
    	return rotmat; 	
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
        
        // TODO should include translation into transformation matrix

        return rotTrans;
    }
    
    public static Matrix4d superposeWithTranslation(Point3d[] x, Point3d[] y) {
        //superpose x onto y
             
        // translate to origin
    	Point3d[] xref = clonePoint3dArray(x);
    	Point3d xtrans = centroid(xref);
    	xtrans.negate();
    	translate(xtrans, xref);

    	Point3d[] yref = clonePoint3dArray(y);
        Point3d ytrans = centroid(yref);
        ytrans.negate();
        translate(ytrans, yref); 
        
        // calculate rotational component (rotation around origin)
        Quat4d q = quaternionOrientation(xref, yref);
        q.conjugate();
        Matrix4d rotTrans = new Matrix4d();
        rotTrans.set(q);   
 
        // combine with x -> origin translation
        Matrix4d trans = new Matrix4d();
        trans.setIdentity();
        trans.setTranslation(new Vector3d(xtrans));
        rotTrans.mul(rotTrans, trans);

        // combine with origin -> y translation
        ytrans.negate();  
        Matrix4d transInverse = new Matrix4d(); 
        transInverse.setIdentity();     
        transInverse.setTranslation(new Vector3d(ytrans));
        rotTrans.mul(transInverse, rotTrans);
        
        // transform x coordinates onto y coordinate frame
        transform(rotTrans, x);

        return rotTrans;
    }
    
    public static Matrix4d superposeAtOrigin(Point3d[] x, Point3d[] y) {
        Quat4d q = quaternionOrientation(x, y);
        q.conjugate();
        
        Matrix4d rotTrans = new Matrix4d();
        rotTrans.set(q);
        
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
        center.scale(1.0/x.length);
        return center;
    }
    
    private static Quat4d quaternionOrientation(Point3d[] a, Point3d[] b)  {
        Matrix m = calcFormMatrix(a, b);
        EigenvalueDecomposition eig = m.eig();
        double[][] v = eig.getV().getArray();
        Quat4d q = new Quat4d(v[1][3], v[2][3], v[3][3], v[0][3]);
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
    public void calcQCPSuperpositionRotationOnly(Point3d[] x, Point3d[] y) {
        double[] a = new double[9];
  //      long t1 = System.nanoTime();
        double E0 = innerProduct(a, y, x);
  //      long t2 = System.nanoTime();
    //    System.out.println("inter product: " + (t2-t1));
        
        /* calculate the RMSD & rotational matrix */
        double threshold = -1;
        if (!rmsdOnly) {
        	threshold = maxRmsd;
        }
        
     //   long t3 = System.nanoTime();
        int ret = fastCalcRMSDAndRotation(a, E0, x.length, threshold);
     //   long t4 = System.nanoTime();
     //   System.out.println("calc rmsd + rot: " + (t4-t3));
     //   long t5 = System.nanoTime();
        transformation.set(rotmat,new Vector3d(0,0,0), 1);
    //    long t6 = System.nanoTime();
     //   System.out.println("create 4x4: " + (t6-t5));
 
    }
  
    /* Superposition coords2 onto coords1 -- in other words, coords2 is rotated, coords1 is held fixed */
    public void calcRMSDRotationalMatrix(Point3d[] x, Point3d[] y) {
        //superpose x onto y
        
        // translate to origin
    	Point3d[] xref = clonePoint3dArray(x);
    	Point3d xtrans = centroid(xref);
  //  	System.out.println("x centroid: " + xtrans);
    	xtrans.negate();
    	translate(xtrans, xref);

    	Point3d[] yref = clonePoint3dArray(y);
        Point3d ytrans = centroid(yref);
   // 	System.out.println("y centroid: " + ytrans);
        ytrans.negate();
        translate(ytrans, yref); 
;

        /* calculate the (weighted) inner product of two structures */
        double[] a = new double[9];
        double E0 = innerProduct(a, yref, xref);
        
        /* calculate the RMSD & rotational matrix */
        double threshold = -1;
        if (!rmsdOnly) {
        	threshold = maxRmsd;
        }
        
        int ret = fastCalcRMSDAndRotation(a, E0, x.length, threshold);
 //       long t1 = System.nanoTime();
        transformation.set(rotmat,new Vector3d(0,0,0), 1);
  //      long t2 = System.nanoTime();
  //      System.out.println("create transformation: " + (t2-t1));
 //       System.out.println("m3d -> m4d");
 //       System.out.println(transformation);

//        // combine with x -> origin translation
//        Matrix4d trans = new Matrix4d();
//        trans.setIdentity();
//        trans.setTranslation(new Vector3d(xtrans));
//        transformation.mul(transformation, trans);
//        System.out.println("setting xtrans");
//        System.out.println(transformation);
//
//        // combine with origin -> y translation
//        ytrans.negate();  
//        Matrix4d transInverse = new Matrix4d(); 
//        transInverse.setIdentity();     
//        transInverse.setTranslation(new Vector3d(ytrans));
//        transformation.mul(transInverse, transformation);
//        System.out.println("setting ytrans");
//        System.out.println(transformation);
        
        // transform x coordinates onto y coordinate frame
  //      transform(transformation, xref);

  //      System.out.println("QCP rmsd: " + rmsd + " explicit rmsd: " + rmsd(xref, y));
  //      for (int i = 0; i < x.length; i++) {
   //     	System.out.println(xref[i] +" - " + y[i]);
    //    }
    }
    
    /** 
     * http://theobald.brandeis.edu/qcp/qcprot.c
     * @param A
     * @param coords1
     * @param coords2
     * @return
     */
    private double innerProduct(double[] A, Point3d[] coords1, Point3d[] coords2) {
    	double          x1, x2, y1, y2, z1, z2;
        int             i;
   //     const double   *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2];
   //     const double   *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2];
        double          G1 = 0.0, G2 = 0.0;

        A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

        if (weight != null)
        {
            for (i = 0; i < coords1.length; ++i)
            {
                x1 = weight[i] * coords1[i].x;
                y1 = weight[i] * coords1[i].y;
                z1 = weight[i] * coords1[i].z;

                G1 += x1 * coords1[i].x + y1 * coords1[i].y + z1 * coords1[i].z;

                x2 = coords2[i].x;
                y2 = coords2[i].y;
                z2 = coords2[i].z;

                G2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

                A[0] +=  (x1 * x2);
                A[1] +=  (x1 * y2);
                A[2] +=  (x1 * z2);

                A[3] +=  (y1 * x2);
                A[4] +=  (y1 * y2);
                A[5] +=  (y1 * z2);

                A[6] +=  (z1 * x2);
                A[7] +=  (z1 * y2);
                A[8] +=  (z1 * z2);   
            }
        }
        else
        {
            for (i = 0; i < coords1.length; ++i)
            {
                x1 = coords1[i].x;
                y1 = coords1[i].y;
                z1 = coords1[i].z;
                G1 += x1 * x1 + y1 * y1 + z1 * z1;

                x2 = coords2[i].x;
                y2 = coords2[i].y;
                z2 = coords2[i].z;

                G2 += (x2 * x2 + y2 * y2 + z2 * z2);

                A[0] +=  (x1 * x2);
                A[1] +=  (x1 * y2);
                A[2] +=  (x1 * z2);

                A[3] +=  (y1 * x2);
                A[4] +=  (y1 * y2);
                A[5] +=  (y1 * z2);

                A[6] +=  (z1 * x2);
                A[7] +=  (z1 * y2);
                A[8] +=  (z1 * z2);  
            }
        }

        return (G1 + G2) * 0.5;
    }

    private int fastCalcRMSDAndRotation(double[] A, double E0, int len, double minScore)
    {
        double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
        double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
               SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
               SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
               SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
        double[] C = new double[4];
        double[] rot = new double[9];
        int i;
        double mxEigenV; 
        double oldg = 0.0;
        double b, a, delta, rms, qsqr;
        double q1, q2, q3, q4, normq;
        double a11, a12, a13, a14, a21, a22, a23, a24;
        double a31, a32, a33, a34, a41, a42, a43, a44;
        double a2, x2, y2, z2; 
        double xy, az, zx, ay, yz, ax; 
        double a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132; 
        double evecprec = 1d-6;
        double evalprec = 1d-11;

 //       long t1 = System.nanoTime();
        Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
        Syx = A[3]; Syy = A[4]; Syz = A[5];
        Szx = A[6]; Szy = A[7]; Szz = A[8];

        Sxx2 = Sxx * Sxx;
        Syy2 = Syy * Syy;
        Szz2 = Szz * Szz;

        Sxy2 = Sxy * Sxy;
        Syz2 = Syz * Syz;
        Sxz2 = Sxz * Sxz;

        Syx2 = Syx * Syx;
        Szy2 = Szy * Szy;
        Szx2 = Szx * Szx;

        SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
        Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

        C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
        C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

        SxzpSzx = Sxz + Szx;
        SyzpSzy = Syz + Szy;
        SxypSyx = Sxy + Syx;
        SyzmSzy = Syz - Szy;
        SxzmSzx = Sxz - Szx;
        SxymSyx = Sxy - Syx;
        SxxpSyy = Sxx + Syy;
        SxxmSyy = Sxx - Syy;
        Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

        C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

        mxEigenV = E0;      
 
        for (i = 0; i < 50; ++i)
        {
            oldg = mxEigenV;
            x2 = mxEigenV*mxEigenV;
            b = (x2 + C[2])*mxEigenV;
            a = b + C[1];
            delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
            mxEigenV -= delta;

            if (Math.abs(mxEigenV - oldg) < Math.abs(evalprec*mxEigenV))
                break;
        }

        if (i == 50) 
           System.err.println("More than %d iterations needed!" + i);

        /* the fabs() is to guard against extremely small, but *negative* numbers due to floating point error */
        rms = Math.sqrt(Math.abs(2.0 * (E0 - mxEigenV)/len));
        rmsd = rms;

  //      long t2 = System.nanoTime();
  //      System.out.println("QCP setup: " + (t2-t1));
        
   //     long t3 = System.nanoTime();
        if (minScore > 0) 
            if (rms < minScore)
                return (-1); // Don't bother with rotation. 

        a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx;
        a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
        a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
        a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
        a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
        a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
        a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
        q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
        q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
        q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
        q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

        qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;
        
 //       long t4 = System.nanoTime();

  //      System.out.println("cal rot: " + (t4-t3));

    /* The following code tries to calculate another column in the adjoint matrix when the norm of the 
       current column is too small.
       Usually this commented block will never be activated.  To be absolutely safe this should be
       uncommented, but it is most likely unnecessary.  
    */
        
 //       long t5 = System.nanoTime();
        if (qsqr < evecprec)
        {
            q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
            q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
            q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
            q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

            if (qsqr < evecprec)
            {
                double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
                double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
                double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

                q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
                q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
                q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
                q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
                qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

                if (qsqr < evecprec)
                {
                    q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                    q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                    q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                    q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                    qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                    
                    if (qsqr < evecprec)
                    {
                        /* if qsqr is still too small, return the identity matrix. */
                        rot[0] = rot[4] = rot[8] = 1.0;
                        rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;

                        return 0;
                    }
                }
            }
        }
 //       long t6 = System.nanoTime();
 //       System.out.println("special case: " + (t6-t5));

 //       long t7 = System.nanoTime();
        normq = Math.sqrt(qsqr);
        q1 /= normq;
        q2 /= normq;
        q3 /= normq;
        q4 /= normq;

        a2 = q1 * q1;
        x2 = q2 * q2;
        y2 = q3 * q3;
        z2 = q4 * q4;

        xy = q2 * q3;
        az = q1 * q4;
        zx = q4 * q2;
        ay = q1 * q3;
        yz = q3 * q4;
        ax = q1 * q2;

        rot[0] = a2 + x2 - y2 - z2;
        rot[1] = 2 * (xy + az);
        rot[2] = 2 * (zx - ay);
        rot[3] = 2 * (xy - az);
        rot[4] = a2 - x2 + y2 - z2;
        rot[5] = 2 * (yz + ax);
        rot[6] = 2 * (zx + ay);
        rot[7] = 2 * (yz - ax);
        rot[8] = a2 - x2 - y2 + z2;

  //      long t8 = System.nanoTime();
  //      System.out.println("rot elements: " + (t8-t7));
        
 //       long t9 = System.nanoTime();
        rotmat.set(rot);
 //       long t10 = System.nanoTime();
  //      System.out.println("matrix3d: " + (t10-t9));

        
    //    System.out.println("New RMSD: " + rmsd);
  //      Matrix3d m = new Matrix3d(rot);
     //   System.out.println("New rotation matrix");
     //   System.out.println(m);
   //     long t11 = System.nanoTime();
   //     System.out.println("total qcp: " + (t11-t1));
        return 1;
    }
	
}
