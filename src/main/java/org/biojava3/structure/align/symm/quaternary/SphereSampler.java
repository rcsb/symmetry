package org.biojava3.structure.align.symm.quaternary;

import java.util.ArrayList;
import java.util.List;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.AxisAngle4f;
import javax.vecmath.Point3i;
import javax.vecmath.Quat4d;
import javax.vecmath.Tuple3i;

// Generate the permutations and sign changes for a Triple.
class Permute {

	private List<Point3i> triples = new ArrayList<Point3i>();

	Permute(Point3i t) {
		// assert(x.a >= x.b && x.b >= x.c && x.c >= 0);
		// m_arr.push_back(x);
		Point3i tmp = new Point3i();
		tmp.x = t.x;
		tmp.y = t.y;
		tmp.z = t.z;
		triples.add(tmp);
		int n = 1;
		// Do the sign changes
		// if (x.a != 0) {
		// for (size_t i = 0; i < n; ++i)
		// m_arr.push_back(Triple(-m_arr[i].a, m_arr[i].b, m_arr[i].c));
		// n *= 2;
		// }
		if (t.x != 0) {
			for (int i = 0; i < n; ++i) {
				Tuple3i m = triples.get(i);
			
				triples.add(new Point3i(-m.x, m.y, m.z));
			}
			n *= 2;
		}
		// if (t.b != 0) {
		// for (size_t i = 0; i < n; ++i)
		// m_arr.push_back(Triple(m_arr[i].a, -m_arr[i].b, m_arr[i].c));
		// n *= 2;
		// }
		if (t.y != 0) {
			for (int i = 0; i < n; ++i) {
				Point3i m = triples.get(i);
				triples.add(new Point3i(m.x, -m.y, m.z));
			}
			n *= 2;
		}
		// if (x.c != 0) {
		// for (size_t i = 0; i < n; ++i)
		// m_arr.push_back(Triple(m_arr[i].a, m_arr[i].b, -m_arr[i].c));
		// n *= 2;
		// }
		if (t.z != 0) {
			for (int i = 0; i < n; ++i) {
				Point3i m = triples.get(i);
				triples.add(new Point3i(m.x, m.y, -m.z));
			}
			n *= 2;
		}
		if (t.x == t.y && t.y == t.z) {
			return;
		}
		// With at least two distinct indices we can rotate the set thru 3
		// permuations.
		// for (size_t i = 0; i < n; ++i) {
		// m_arr.push_back(Triple(m_arr[i].b, m_arr[i].c, m_arr[i].a));
		// m_arr.push_back(Triple(m_arr[i].c, m_arr[i].a, m_arr[i].b));
		// }
		// n *= 3;
		for (int i = 0; i < n; ++i) {
			Point3i m = triples.get(i);
			triples.add(new Point3i(m.y, m.z, m.x));
			triples.add(new Point3i(m.z, m.x, m.y));
		}
		n *= 3;

		if (t.x == t.y || t.y == t.z) {
			return;
		}
		// With three distinct indices we can in addition interchange the
		// first two indices (to yield all 6 permutations of 3 indices).
		// for (size_t i = 0; i < n; ++i) {
		// m_arr.push_back(Triple(m_arr[i].b, m_arr[i].a, m_arr[i].c));
		// }
		// n *= 2;
		for (int i = 0; i < n; ++i) {
			Point3i m = triples.get(i);
			triples.add(new Point3i(m.y, m.x, m.z));
		}
		n *= 2;
	}

	public int size() {
		return triples.size();
	}

	public Point3i get(int i) {
		return triples.get(i);
	}
};

/**
 * 
 * @author Peter
 */
public final class SphereSampler {
	private static List<Quat4d> orientations = new ArrayList<Quat4d>();

	// The rotational symmetries of the cube. (Not normalized, since
	// PackSet.Add does this.)
	private static double cubeSyms[][] = {
			{ 1, 0, 0, 0 },
			// 180 deg rotations about 3 axes
			{ 0, 1, 0, 0 },
			{ 0, 0, 1, 0 },
			{ 0, 0, 0, 1 },
			// +/- 120 degree rotations about 4 leading diagonals
			{ 1, 1, 1, 1 }, { 1, 1, 1, -1 }, { 1, 1, -1, 1 }, { 1, 1, -1, -1 },
			{ 1, -1, 1, 1 }, { 1, -1, 1, -1 },
			{ 1, -1, -1, 1 },
			{ 1, -1, -1, -1 },
			// +/- 90 degree rotations about 3 axes
			{ 1, 1, 0, 0 }, { 1, -1, 0, 0 }, { 1, 0, 1, 0 }, { 1, 0, -1, 0 },
			{ 1, 0, 0, 1 }, { 1, 0, 0, -1 },
			// 180 degree rotations about 6 face diagonals
			{ 0, 1, 1, 0 }, { 0, 1, -1, 0 }, { 0, 1, 0, 1 }, { 0, 1, 0, -1 },
			{ 0, 0, 1, 1 }, { 0, 0, 1, -1 }, };
//	 private static double delta = 0.25970;
//	 private static double sigma = 0.00;
//	 private static int ntot = 1992;
//	 private static int ncell = 83;
//	 private static int nent = 7;
//	 private static double maxrad = 16.29;
//	 private static double coverage = 2.42065;
//	 private static int[] k = {0, 1, 2, 2, 2, 3, 3};
//	 private static int[] l = {0, 1, 0, 2, 2, 1, 3};
//	 private static int[] m = {0, 1, 0, 0, 2, 1, 1};
//	 private static double[] w = {1.665264, 1.517726, 1.489794, 1.205193,
//	 1.146566, 0.973349, 0.552456};
//	 private static int[] mult = {1, 8, 6, 12, 8, 24, 24};

	// # Orientation set c48u83, number = 1992, radius = 16.29 degrees
	// # $Id: c48u83.grid 6102 2006-02-21 19:45:40Z ckarney $
	// # For more information, See http://charles.karney.info/orientation/
	// format grid
	// 0.25970 0.00 1992 83 7 16.29 2.42065
	// 0 0 0 1.665264 16.29 1
	// 1 1 1 1.517726 16.29 8
	// 2 0 0 1.489794 16.29 6
	// 2 2 0 1.205193 16.00 12
	// 2 2 2 1.146566 15.33 8
	// 3 1 1 0.973349 16.29 24
	// 3 3 1 0.552456 14.88 24

//	private static double delta = 0.19415;
//	private static double sigma = 0.00;
//	private static int ntot = 4344;
//	private static int ncell = 181;
//	private static int nent = 13;
//	private static double maxrad = 12.29;
//	private static double coverage = 2.27013;
//	private static int[] k = { 0, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4 };
//	private static int[] l = { 0, 1, 0, 2, 2, 1, 3, 3, 0, 2, 2, 4, 4 };
//	private static int[] m = { 0, 1, 0, 0, 2, 1, 1, 3, 0, 0, 2, 0, 2 };
//	private static double[] w = { 1.557196, 1.475638, 1.449892, 1.353190,
//			1.265726, 1.257478, 1.143450, 0.867360, 0.916638, 0.941799,
//			0.876165, 0.680232, 0.446640 };
//	private static int[] mult = { 1, 8, 6, 12, 8, 24, 24, 8, 6, 24, 24, 12, 24 };
	// # Orientation set c48u181, number = 4344, radius = 12.29 degrees
	// # $Id: c48u181.grid 6102 2006-02-21 19:45:40Z ckarney $
	// # For more information, See http://charles.karney.info/orientation/
	// format grid
	// 0.19415 0.00 4344 181 13 12.29 2.27013
	// 0 0 0 1.557196 12.29 1
	// 1 1 1 1.475638 12.29 8
	// 2 0 0 1.449892 12.29 6
	// 2 2 0 1.353190 12.16 12
	// 2 2 2 1.265726 11.86 8
	// 3 1 1 1.257478 11.99 24
	// 3 3 1 1.143450 12.29 24
	// 3 3 3 0.867360 11.12 8
	// 4 0 0 0.916638 11.70 6
	// 4 2 0 0.941799 11.61 24
	// 4 2 2 0.876165 12.29 24
	// 4 4 0 0.680232 10.95 12
	// 4 4 2 0.446640 12.25 24
	
	private static double delta = 0.15846;
	private static double sigma = 0.00;
	private static int ntot = 7416;
	private static int ncell = 309;
	private static int nent = 18;
//	private static double maxrad = 10.07;
//	private static double coverage = 2.13338;
	private static int[] k = { 0, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5};
	private static int[] l = { 0, 1, 0, 2, 2, 1, 3, 3, 0, 2, 2, 4, 4, 4, 1, 3, 3, 5};
	private static int[] m = { 0, 1, 0, 0, 2, 1, 1, 3, 0, 0, 2, 0, 2, 4, 1, 1, 3, 1};
//	private static double[] w = { 1.557196, 1.475638, 1.449892, 1.353190,
//			1.265726, 1.257478, 1.143450, 0.867360, 0.916638, 0.941799,
//			0.876165, 0.680232, 0.446640, 0.0, 0.0, 0.0, 0.0, 0.0};
	private static int[] mult = { 1, 8, 6, 12, 8, 24, 24, 8, 6, 24, 24, 12, 24, 8, 24, 48, 24, 24};
	
//	# Orientation set c48u309, number = 7416, radius = 10.07 degrees
//	# $Id: c48u309.grid 6102 2006-02-21 19:45:40Z ckarney $
//	# For more information, See http://charles.karney.info/orientation/
//	format grid
//	0.15846 0.00   7416  309  18 10.07  2.13338
//	 0  0  0  1.461532  10.07   1
//	 1  1  1  1.409259  10.07   8
//	 2  0  0  1.392463  10.07   6
//	 2  2  0  1.328138  10.00  12
//	 2  2  2  1.268130   9.83   8
//	 3  1  1  1.282781   9.90  24
//	 3  3  1  1.172478   9.70  24
//	 3  3  3  1.075693   9.39   8
//	 4  0  0  1.233711  10.07   6
//	 4  2  0  1.121495   9.68  24
//	 4  2  2  1.095831   9.53  24
//	 4  4  0  0.918990   9.38  12
//	 4  4  2  0.969269   9.49  24
//	 4  4  4  0.691509   8.83   8
//	 5  1  1  0.872382  10.07  24
//	 5  3  1  0.757191   9.50  48
//	 5  3  3  0.740975   9.00  24
//	 5  5  1  0.782872   9.50  24


	static {
		createSphereSet();
	}

	// this class cannot be instantiated
	private SphereSampler() {
	};

	public static int getSphereCount() {
		return orientations.size();
	}

	public static Quat4d getQuat4d(int index) {
		return orientations.get(index);
	}

	public static void getAxisAngle(int index, AxisAngle4f axisAngle) {
		axisAngle.set(orientations.get(index));
	}

	public static void getAxisAngle(int index, AxisAngle4d axisAngle) {
		axisAngle.set(orientations.get(index));
	}

	// Convert from index to position. The sinh scaling tries to compensate
	// for the bunching up that occurs when [1 x y z] is projected onto the
	// unit sphere.
	private static double pind(double ind, double delta, double sigma) {
		return (sigma == 0) ? ind * delta : Math.sinh(sigma * ind * delta)
				/ sigma;
	}

	private static void createSphereSet() {
		for (int i = 0; i < IcosahedralSampler.getSphereCount(); i++) {
			orientations.add(IcosahedralSampler.getQuat4d(i));
		}
		int ncell1 = 0;
		for (int n = 0; n < nent; ++n) {
			Permute p = new Permute(new Point3i(k[n], l[n], m[n]));
//			System.out.println("Permutate = " + k[n] + "," + l[n] + "," + m[n]
//					+ ": " + p.size());
			assert (mult[n] == p.size());
			for (int i = 0; i < mult[n]; ++i) {
				Point3i t = p.get(i);
				orientations.add(new Quat4d(1.0, pind(0.5 * t.x, delta, sigma), pind(
						0.5 * t.y, delta, sigma), pind(0.5 * t.z, delta, sigma)));
			}
			ncell1 += mult[n];
		}
		assert (ncell1 == ncell);
//			System.out.println("not fine");
			int nc = orientations.size();
			assert (nc == ncell);
			for (int n = 1; n < 24; ++n) {
				Quat4d q = new Quat4d(cubeSyms[n][0], cubeSyms[n][1],
						cubeSyms[n][2], cubeSyms[n][3]);
//				System.out.println("q: " + q);
				for (int i = 0; i < nc; ++i) {
					 Quat4d qs = new Quat4d();
					 qs.mul(q, orientations.get(i));
					 orientations.add(qs);
	//				s.add(times(q, s.getOrientation(i)), s.getWeight(i));
				}
			}
			assert (orientations.size() == ntot);
			// s.clear();
		// for (int i = 0; i < s.size(); i++) {
		// System.out.println(s.getOrientation(i) + ", " + s.getWeight(i));
		// }
	}
}
