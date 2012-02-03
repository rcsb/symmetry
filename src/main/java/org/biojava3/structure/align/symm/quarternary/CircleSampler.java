package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.List;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.AxisAngle4f;
import javax.vecmath.Point3i;
import javax.vecmath.Quat4d;
import javax.vecmath.Tuple3i;
import javax.vecmath.Vector3d;


/**
 * 
 * @author Peter
 */
public final class CircleSampler {
	private static List<AxisAngle4d> orientations = new ArrayList<AxisAngle4d>();

	static {
		createSphereSet();
	}

	// this class cannot be instantiated
	private CircleSampler() {
	};

	public static int getSphereCount() {
		return orientations.size();
	}

//	public static Quat4d getQuat4d(int index) {
//		return new Quat4d();
//	}

	public static void getAxisAngle(int index, AxisAngle4f axisAngle) {
	    orientations.get(index);
	}

	public static void getAxisAngle(int index, AxisAngle4d axisAngle) {
		axisAngle.set(orientations.get(index));
	}

	private static void createSphereSet() {		
		for (int i = 0; i < 1800; i++) {
			double angle = Math.toRadians(0.1*i);
			Vector3d axis = new Vector3d(0.0, Math.sin(angle), Math.cos(angle));
			AxisAngle4d axisAngle = new AxisAngle4d();
			axisAngle.set(axis, 0.0);
//			System.out.println("circle: " + axis + "axisangle: " + axisAngle);
			orientations.add(axisAngle);
		}	
	}

}
