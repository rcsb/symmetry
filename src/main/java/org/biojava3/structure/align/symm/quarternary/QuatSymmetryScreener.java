package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.vecmath.Matrix4d;
import javax.vecmath.Matrix4f;
import javax.vecmath.Point3d;
import javax.vecmath.Point3f;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;

public class QuatSymmetryScreener {
	private Subunits subunits = null;
	private final double SYMMETRY_THRESHOLD = 0.2;
//	private float distanceThreshold = 2.0f;
	private double distanceThreshold = 3.0; // problem with 2BTV
//	private float distanceThreshold = 4.0; // problem with 1Z8Y
//	private float distanceThreshold = 7.5; // this works with 1Z8Y
	private List<List<Integer>> orbits = new ArrayList<List<Integer>>();
	private int rises;
	private boolean helical = false;
	private double correlation = 0.0;

	public QuatSymmetryScreener(Subunits subunits) {
		this.subunits = subunits;
	}

	public List<List<Integer>> getOrbits() {
		run();
		return orbits;
	}

	public boolean isHelical() {
		run();
		return helical;
	}

	public int getRises() {
		run();
		return rises;
	}

	public double getCorrelation() {
		run();
		return correlation;
	}

	private void run() {
		if (subunits.getSubunitCount() == 0) {
			return;
		}
		if (orbits.size() == 0) {
			calcOrbits();
			helical = isHelical();
		}
	}
	private void calcOrbits() {
		// check if centers are equidistant from the center of mass
		List<Point3d> centers = subunits.getCenters();
		Point3d centroid = new Point3d(); // centroid is origin
		int n = centers.size();
		double[] distance = new double[n];

		for (int i = 0; i < n; i++) {
			distance[i] = centroid.distance(centers.get(i));
		}
		Arrays.sort(distance);
//		System.out.println("Dist: " + Arrays.toString(distance));

		List<Integer> orbit = new ArrayList<Integer>();
		
		double referenceDistance = distance[0];

		for (int i = 0; i < n - 1; i++) {
			int j = i + 1;
			double ratio = distance[j] / referenceDistance; // could use dist Sq -> give ratio Sq
//			System.out.println("i: " + i + " d1 " + distance[j] + "/"
//					+ referenceDistance + " r: " + ratio);
			if (ratio < 1.0 + SYMMETRY_THRESHOLD
					&& distance[j] - referenceDistance < distanceThreshold) {
				orbit.add(i);
			} else {
				orbit.add(i);
				orbits.add(orbit);
				orbit = new ArrayList<Integer>();
				referenceDistance = distance[j];
			}
		}
		orbit.add(n - 1);
		orbits.add(orbit);

		calcHelical(distance);
//		for (List<Integer> orb: orbits) {
//			System.out.println("Orbit: " + orb);
//		}
	}
	
	public int getMaximumSymmetryNumber() {
		run();
		if (orbits == null) {
			return -1;
		}
		int symmetryNumber = Integer.MAX_VALUE;

		for (List<Integer> orbit: orbits) {
			if (orbit.size() > 1) {
				symmetryNumber = Math.min(symmetryNumber, orbit.size());
			}
		}

		if (symmetryNumber == Integer.MAX_VALUE) {
			symmetryNumber = 1;
		}

		return symmetryNumber;
	}
	
	private void calcHelical(double[]  distance) {
		helical = false;
		rises = 0;
		correlation = 0;
		
		// try using only an odd number of distances, so that the distance to the
		// closest subunit equal the radius of the helix
		// ...
		
		// radius squared of helix
		double radiusSq = distance[0]*distance[0];
		
		// max distance squared for subunit farthest from center of mass
		double rMaxSq = distance[distance.length-1]*distance[distance.length-1];
		
		// total height of helix
		double t = 2* Math.sqrt(rMaxSq - radiusSq);
		
//		System.out.println("Helix height: " + t);
			
		// height per subunit
		double rise = t/distance.length;
		// distance threshold for identifying repeating subunit distance
		double threshold = 0.5 * rise;
		
//		System.out.println("helical threshold: " + threshold);
		
		// height of each subunit measured from the subunit closest to the
		// center of mass
		double[] height = new double[distance.length];
		
	//	System.out.println("Unique Heights");
		int count = 0;
		for (int i = 0; i < distance.length; i++) {
			double h = Math.sqrt(distance[i]*distance[i] - radiusSq);
			if (i == 0 || h - height[count-1] > threshold) {
				height[count] = h;
//				System.out.println(h);
				count++;
			} 
		}
		if (count < 2) {
			helical = false;
			rises = count;
			return;
		}
		
		// remove first and last entry, these are often outliers
		height = Arrays.copyOfRange(height, 1, count-1);
		count -=2;
		
		double[] index = new double[count];
		for (int i = 0; i < count; i++) {
		    index[i] = (double)i;
		}

		if (count > 2) {
		   correlation = calcPearson(index, height) ;
		} 
//		System.out.println("#rises: " + count + " r: " + correlation);
		helical = correlation > 0.98;
		rises = count+2;
		
		// exclude short "helices" since they are mostly false positives, and those with 5 or fewer subunits
		if (t < 30 || subunits.getSubunitCount() < 6) {
			helical = false;
		}
	}
	
	/**
	 * Calculate Pearson correlation coefficient
	 * 
	 * @return Pearson correlation coefficient
	 */
	public static double calcPearson(double[] x, double[] y) {
		double sumx = 0;
		double sumy = 0;
		double sumxSq = 0;
		double sumySq = 0;
		double pSum = 0;
		int n = x.length;
		for (int i = 0; i < n; i++) {
			sumx += x[i];
			sumy += y[i];
			sumxSq += x[i] * x[i];
			sumySq += y[i] * y[i];
			pSum = pSum + x[i] * y[i];
		}
		double numerator = pSum - sumx * sumy / n;
		double denominator = Math.sqrt((sumxSq - sumx * sumx / n)
				* (sumySq - sumy * sumy / n));
		return numerator / denominator;
	}



}
