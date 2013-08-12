package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.quaternary.geometry.SuperPosition;

public class HelixSolver {
	private Subunits subunits = null;
	private int fold = 1;
	private HelixLayers helixLayers = new HelixLayers(); 
	private QuatSymmetryParameters parameters = null;
	boolean modified = true;

	public HelixSolver(Subunits subunits, int fold, QuatSymmetryParameters parameters) {
		this.subunits = subunits;
		this.fold = fold;
		this.parameters = parameters;
	}
	
	public HelixLayers getSymmetryOperations() {
		if (modified) {
			solve();
			modified = false;
		}
		return helixLayers;
	}

	private void solve() {
		if (! preCheck()) {
			return;
		}
		
		HelicalRepeatUnit unit = new HelicalRepeatUnit(subunits);
		List<Point3d> repeatUnitCenters = unit.getRepeatUnitCenters();
		List<Point3d[]> repeatUnits = unit.getRepeatUnits();
		Set<List<Integer>> permutations = new HashSet<List<Integer>>();
		Set<Integer> subunitsUsed = new HashSet<Integer>();
		int maxCount = 0;
		double minRise = parameters.getMinimumHelixRise() * fold;  // for n-start helix, the rise must be steeper
		double maxRise = 0;
		Map<Integer[], Integer> interactionMap = unit.getInteractingRepeatUnits();

//		for (Integer[] pair: unit.getInteractingRepeatUnits().keySet()) {
		for (Entry<Integer[], Integer> entry : interactionMap.entrySet()) {
			Integer[] pair = entry.getKey();
			int contacts = entry.getValue();
			Point3d[] h1 = null;
			Point3d[] h2 = null;

			//	System.out.println("************** Comparing: " + pair[0] + "-" + pair[1]);
			// trial superposition of repeat unit pairs to get a seed permutation
			h1 = SuperPosition.clonePoint3dArray(repeatUnits.get(pair[0]));
			h2 = SuperPosition.clonePoint3dArray(repeatUnits.get(pair[1]));
			Matrix4d transformation = SuperPosition.superposeWithTranslation(h1, h2);

			double rmsd = SuperPosition.rmsd(h1, h2);		
			double rise = getRise(transformation, repeatUnitCenters.get(pair[0]), repeatUnitCenters.get(pair[1]));
			double angle = getAngle(transformation);

//			if (parameters.isVerbose()) {
//				System.out.println();
//				System.out.println("Original rmsd: " + rmsd);
//				System.out.println("Original rise: " + rise);
//				System.out.println("Original angle: " + Math.toDegrees(angle));
//			}

			if (rmsd > parameters.getRmsdThreshold()) {
				continue;
			}

			if (Math.abs(rise) < minRise) {
				continue;
			}

			List<Integer> permutation = getPermutation(transformation);
//			if (parameters.isVerbose()) {
//				System.out.println("Permutation: " + permutation);
//			}

			// check permutations

			// don't save redundant permutations
			if (permutations.contains(permutation)) {
				continue;
			}
			permutations.add(permutation);

			// in a helix a repeat unit cannot map onto itself
			boolean valid = true;
			for (int j = 0; j < permutation.size(); j++) {
				if (permutation.get(j) == j) {
					valid = false;
					break;
				}
			}
			if (!valid) {
				continue;
			}

			// keep track of which subunits are permutated
			Set<Integer> permSet = new HashSet<Integer>();
			int count = 0;
			for (int i = 0; i < permutation.size(); i++) {
				if (permutation.get(i) != -1) {
					permSet.add(permutation.get(i));
					permSet.add(i);
					count++;
				}
			}
			if (count == 1) {
				continue;
			}

			// superpose all permuted subunits
			List<Point3d> point1 = new ArrayList<Point3d>();
			List<Point3d> point2 = new ArrayList<Point3d>();
			List<Point3d> centers = subunits.getOriginalCenters();
			for (int j = 0; j < permutation.size(); j++) {
				if (permutation.get(j) == -1) {
					continue;
				}
				point1.add(new Point3d(centers.get(j)));
				point2.add(new Point3d(centers.get(permutation.get(j))));
			}

			// for helical symmetry, the number of permuted subunits must be less
			// than the total number of subunits, if they are equal cyclic symmetry was found
			if (point1.size() == subunits.getSubunitCount()) {
				continue;
			}

			h1 = new Point3d[point1.size()];
			h2 = new Point3d[point2.size()];
			point1.toArray(h1);
			point2.toArray(h2);

			// calculate subunit rmsd if at least 3 subunits are available
			double subunitRmsd = 0;
			if (point1.size() > 2) {
				transformation = SuperPosition.superposeWithTranslation(h1, h2);

				subunitRmsd = SuperPosition.rmsd(h1, h2);
				rise = getRise(transformation, repeatUnitCenters.get(pair[0]), repeatUnitCenters.get(pair[1]));
				angle = getAngle(transformation);

				if (parameters.isVerbose()) {
					System.out.println("Subunit rmsd: " + subunitRmsd);
					System.out.println("Subunit rise: " + rise);
					System.out.println("Subunit angle: " + Math.toDegrees(angle));
				}

				if (subunitRmsd > parameters.getRmsdThreshold()) {
					continue;
				}

				if (Math.abs(rise) < minRise) {
					continue;
				}
			}

			// superpose all C alpha traces
			point1.clear();
			point2.clear();
			List<Point3d[]> traces = subunits.getTraces();
			for (int j = 0; j < permutation.size(); j++) {
				if (permutation.get(j) == -1) {
					continue;
				}
				for (Point3d p: traces.get(j)) {
					point1.add(new Point3d(p));
				}
				for (Point3d p: traces.get(permutation.get(j))) {
					point2.add(new Point3d(p));
				}
			}

			h1 = new Point3d[point1.size()];
			h2 = new Point3d[point2.size()];
			point1.toArray(h1);
			point2.toArray(h2);

			transformation = SuperPosition.superposeWithTranslation(h1, h2);

			double traceRmsd = SuperPosition.rmsd(h1, h2);
			
			// TM score should be calcualted separately for each subunit!
			double traceTmScore = SuperPosition.TMScore(h1, h2, h1.length);
			rise = getRise(transformation, repeatUnitCenters.get(pair[0]), repeatUnitCenters.get(pair[1]));
			angle = getAngle(transformation);

//			if (parameters.isVerbose()) {
//				System.out.println("Trace rmsd: " + traceRmsd);
//				System.out.println("Trace rise: " + rise);
//				System.out.println("Trace angle: " + Math.toDegrees(angle));
//			}

			if (traceRmsd > parameters.getRmsdThreshold()) {
				continue;
			}

			if (Math.abs(rise) < minRise) {
				continue;
			}

			// This prevents translational repeats to be counted as helices
			// TODO fine-tune this threshold and add to QuatSymmeteryParameters
			if (angle < Math.toRadians(parameters.getMinimumHelixAngle())) {
				continue;
			}

			// save this helix rot-tranlation
			Helix helix = new Helix();
			helix.setTransformation(transformation);
			helix.setPermutation(permutation);
			helix.setRise(rise);
			helix.setSubunitRmsd(subunitRmsd);
			helix.setTraceRmsd(traceRmsd);
			helix.setTraceTmScoreMin(traceTmScore);
			helix.setFold(fold);
			helix.setContacts(contacts);
			helix.setRepeatUnits(unit.getRepeatUnitIndices());
			helixLayers.addHelix(helix);

			subunitsUsed.addAll(permSet);
			maxCount = Math.max(maxCount, count);
			maxRise = Math.max(maxRise, rise);
		}

		// if too few subunits are permuted, it's likely an invalid solution
		if (maxCount < (subunits.getSubunitCount()-this.fold)*0.8) {
//			if (parameters.isVerbose()) {
//				System.out.println("maxCount too low: " + maxCount + " out of " + subunits.getSubunitCount());
//			}
			helixLayers = new HelixLayers();
		}

		// if not all subunits have been used, there can't be full helical symmetry
		if (subunitsUsed.size() != subunits.getSubunitCount()) {
			helixLayers = new HelixLayers();
		}

		return;
	}

	private boolean preCheck() {
		if (subunits.getSubunitCount() < 3) {
			return false;
		}
		List<Integer> folds = this.subunits.getFolds();
		int maxFold = folds.get(folds.size()-1);
		return maxFold > 1;
	}

	/** 
	 * Returns a permutation of subunit indices for the given helix transformation.
	 * An index of -1 is used to indicate subunits that do not superpose onto any other subunit.
	 * @param transformation
	 * @return
	 */
	private List<Integer> getPermutation(Matrix4d transformation) {
		double rmsdThresholdSq = Math.pow(this.parameters.getRmsdThreshold(), 2);

		List<Point3d> centers = subunits.getOriginalCenters();

		List<Integer> permutations = new ArrayList<Integer>(centers.size());
		double[] dSqs = new double[centers.size()];
		boolean[] used = new boolean[centers.size()];
		Arrays.fill(used, false);

		for (Point3d center: centers) {
			Point3d tCenter = new Point3d(center);
			transformation.transform(tCenter);
			int permutation = -1;
			double minDistSq = Double.MAX_VALUE;
			for (int i = 0; i < centers.size(); i++) {
				if (! used[i]) {
					double dSq = tCenter.distanceSquared(centers.get(i));
					if (dSq < minDistSq && dSq <= rmsdThresholdSq) {
						minDistSq = dSq;
						permutation = i; 
						dSqs[i] = dSq;
					}
				}
			}
			// can't map to itself
			if (permutations.size() == permutation) {
				permutation = -1;
			}

			if (permutation != -1) {
				used[permutation] = true;
			}

			permutations.add(permutation);
		}

		return permutations;
	}

	/**
	 * Returns the rise of a helix given the subunit centers of two adjacent
	 * subunits and the helix transformation
	 * @param transformation helix transformation
	 * @param p1 center of one subunit
	 * @param p2 center of an adjacent subunit
	 * @return
	 */
	private static double getRise(Matrix4d transformation, Point3d p1, Point3d p2) {
		AxisAngle4d axis = getAxisAngle(transformation);
		Vector3d h = new Vector3d(axis.x, axis.y, axis.z);
		Vector3d p = new Vector3d();
		p.sub(p1, p2);
		return p.dot(h);
	}
	
	/**
	 * Returns the pitch angle of the helix
	 * @param transformation helix transformation
	 * @return
	 */
	private static double getAngle(Matrix4d transformation) {
		return getAxisAngle(transformation).angle;
	}
	
	/**
	 * Returns the AxisAngle of the helix transformation
	 * @param transformation helix transformation
	 * @return
	 */
	private static AxisAngle4d getAxisAngle(Matrix4d transformation) {
		AxisAngle4d axis = new AxisAngle4d();
		axis.set(transformation);
		return axis;
	}
}
