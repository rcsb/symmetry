package org.biojava3.structure.quaternary.misc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.quaternary.core.ChainClusterer;
import org.biojava3.structure.quaternary.core.FindQuarternarySymmetry;
import org.biojava3.structure.quaternary.core.PermutationGenerator;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.core.Subunits;
import org.biojava3.structure.quaternary.geometry.SuperPosition;

public class SubunitMapper {
	private QuatSymmetryParameters parameters = new QuatSymmetryParameters();
	private Structure structure1 = null;
	private Structure structure2 = null;
	
	private RotationGroup rotationGroup1 = null;
	private RotationGroup rotationGroup2 = null;
	private Subunits subunits1 = null; 
	private Subunits subunits2 = null;
	
	public SubunitMapper(Structure structure1, Structure structure2) {
		this.structure1 = structure1;
		this.structure2 = structure2;
	}
	
	public List<List<Integer>> getMappings() {
		List<List<Integer>> mappings = new ArrayList<List<Integer>>();
		if (!setupSymmetryInfo()) {
			return mappings;
		}
		
		if (subunits1.getCenters().size() <= 4) {
			mappings = getAllMappings();
			if (! mappings.isEmpty()) {
				return mappings;
			}
		}
		
		long t1 = System.nanoTime();
		List<Integer> map1 = getMappingFromAxisAlignment(1);
		long t2 = System.nanoTime();
		System.out.println("Axix mapping 1: " + (t2-t1)/1000000 + " ms");
		System.out.println("____________________Axis mapping 2 _____________________________");
		List<Integer> map2 = getMappingFromAxisAlignment(-1);
		long t3 = System.nanoTime();
		System.out.println("Axix mapping 2: " + (t3-t2)/1000000 + " ms");
		if (! map1.isEmpty()) {
			mappings.add(map1);
		}
		if (! map2.isEmpty()) {
			mappings.add(map2);
		}
		List<Integer> map3 = getMappingByStructuralSuperposition();
		if (! map3.isEmpty()) {
			mappings.add(map3);
		}
		return mappings;
	}
	
	private boolean setupSymmetryInfo() {
		FindQuarternarySymmetry finder1 = getQuaternarySymmetry(structure1);
		rotationGroup1 = finder1.getRotationGroup();	
		subunits1 = finder1.getSubunits();
		finder1 = null; // gc
		int n1 = subunits1.getCenters().size();	
		String p1 = rotationGroup1.getPointGroup();

		FindQuarternarySymmetry finder2 = getQuaternarySymmetry(structure2);
		rotationGroup2 = finder2.getRotationGroup();	
		subunits2 = finder2.getSubunits();
		finder2 = null; // gc
		int n2 = subunits2.getCenters().size();
		String p2 = rotationGroup2.getPointGroup();
		
		if (n1 != n2) {
			System.out.println("Cannot superpose structures with different number of subunits: " + n1 + " vs. " + n2);
			return false;
		}
		if (!p1.equals(p2)) {
			System.out.println("Cannot superpose structures with different point groups: " + p1 + " vs. " + 22);
			return false;
		}
		return true;
	}
	
	public List<List<Integer>> getAllMappings() {
		List<List<Integer>> mappings = new ArrayList<List<Integer>>();
		PermutationGenerator g = new PermutationGenerator(subunits1.getCenters().size());
		while (g.hasMore()) {
			//         for (List<Integer> mapping: mappings) {
			int[] p = g.getNext();
			List<Integer> mapping = new ArrayList<Integer>();
			for (int e: p) {
				mapping.add(e);
			}
			mappings.add(mapping);
		}
		return mappings;
	}
	
	/**
	 * Calculate approximate relative orientation of subunits
	 * in the two structures by superposing the principal axis
	 * and one reference subunit
	 * @return
	 */
	private List<Integer> getMappingFromAxisAlignment(int sign) {
		Point3d[] points1 = getAxisPoints(subunits1, rotationGroup1, 1);
		Point3d[] points2 = getAxisPoints(subunits2, rotationGroup2, sign);

		AxisAngle4d axis = new AxisAngle4d();
		Matrix4d transformation = SuperPosition.superposeAtOrigin(points1, points2, axis);
		System.out.println("axis superposition:  rmsd: " + SuperPosition.rmsd(points1, points2));
		System.out.println(Arrays.toString(points1));
		System.out.println(Arrays.toString(points2));
		
		List<Vector3d> vectors = subunits1.getUnitVectors();
		List<Point3d> transformed = new ArrayList<Point3d>();
		for (Vector3d v: vectors) {
			Point3d pt = new Point3d(v);
			transformation.transform(pt);
			transformed.add(pt);
		}
		
		List<Point3d> original = new ArrayList<Point3d>();
		for (Vector3d v: subunits2.getUnitVectors()) {
		      original.add(new Point3d(v));
		}
		
		System.out.println("Superposed s1/s2");
		System.out.println(transformed);
		System.out.println(original);
		
		List<Integer> mapping = calcMapping(original, transformed);
		
		return mapping;
	}
	
	private List<Integer> getMappingByStructuralSuperposition() {
		StructureAlignment algo = null;
		try {
			algo = StructureAlignmentFactory.getAlgorithm(SmithWaterman3Daligner.algorithmName);
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ChainClusterer seqGroup1 = new ChainClusterer(structure1, parameters);
		List<Atom[]> cas1 = seqGroup1.getCalphaTraces();
		ChainClusterer seqGroup2 = new ChainClusterer(structure2, parameters);
		List<Atom[]> cas2 = seqGroup2.getCalphaTraces();
		
	
//		System.out.println("Chain1: orig:");
//		for (Atom a: cas1.get(0)) {
//			System.out.println(Arrays.toString(a.getCoords()));
//		}
//		
//		System.out.println("Chain2: orig:");
//		for (Atom a: cas2.get(0)) {
//			System.out.println(Arrays.toString(a.getCoords()));
//		}
		
		AFPChain afpChain = null;
		try {
			afpChain = algo.align(cas2.get(0), cas1.get(0));
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("Identity: " + afpChain.getIdentity());
		System.out.println("Chain rmsd: " + afpChain.getChainRmsd());
		
		Atom[] trans = afpChain.getBlockShiftVector();
		Matrix[] mb = afpChain.getBlockRotationMatrix();	
		
		if (mb == null || mb[0] == null || trans == null | trans[0] == null) {
			return new ArrayList<Integer>();
		}
		
		
		List<Point3d> transformed = new ArrayList<Point3d>();
		for (int i = 0; i < subunits1.getOriginalCenters().size(); i++) {
			Point3d p = new Point3d(subunits1.getOriginalCenters().get(i));
			double[] c = {p.x, p.y, p.z};
			Atom a = new AtomImpl();
			a.setCoords(c);
			Calc.rotate(a, mb[0]);
			Calc.shift(a, trans[0]);
			c = a.getCoords();
			p.set(c);
			transformed.add(p);
		}
		
		List<Integer> map = calcMapping(subunits2.getOriginalCenters(), transformed);

		return map;
	}
	
	private static Vector3d[] getAxisVectors(Subunits subunits, RotationGroup rotGroup, int sign) {
		Vector3d[] vectors = new Vector3d[3];
		// add points on principal axis 10 A in each direction
		AxisAngle4d axisAngle = rotGroup.getRotation(0).getAxisAngle();
		vectors[0] = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		vectors[0].normalize();
		vectors[0].scale(sign);
		vectors[1] = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		vectors[1].normalize();
		vectors[1].scale(-sign);
		
		// find other axis that is perpendicular to principal axis
		for (int i = 0; i < rotGroup.getOrder(); i++) {
			System.out.println("dir: " + rotGroup.getRotation(i).getDirection() + " fold: " + rotGroup.getRotation(i).getFold() + " angle: " + rotGroup.getRotation(i).getAxisAngle().angle);
			if (rotGroup.getRotation(i).getDirection() == 1) {
				axisAngle = rotGroup.getRotation(i).getAxisAngle();
				vectors[2] = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				vectors[2].normalize();
				vectors[2].scale(1);
				break;
			}
		}
		
//		vectors[2] = null;
		
		if (vectors[2] == null) {
			Vector3d v =  new Vector3d(subunits.getCenters().get(0));
			v.normalize();
			v.scale(1);
            vectors[2] = v;
		}
		return vectors;
	}
	
	private static Point3d[] getAxisPoints(Subunits subunits, RotationGroup rotGroup, int sign) {
		Point3d[] points = new Point3d[3];
		// add points on principal axis 10 A in each direction
		AxisAngle4d axisAngle = rotGroup.getRotation(0).getAxisAngle();
		points[0] = new Point3d(axisAngle.x, axisAngle.y, axisAngle.z);
		points[0].scale(sign);
		points[1] = new Point3d(axisAngle.x, axisAngle.y, axisAngle.z);
		points[1].scale(-sign);
		
		// find other axis that is perpendicular to principal axis
		for (int i = 0; i < rotGroup.getOrder(); i++) {
			System.out.println("dir: " + rotGroup.getRotation(i).getDirection() + " fold: " + rotGroup.getRotation(i).getFold() + " angle: " + rotGroup.getRotation(i).getAxisAngle().angle);
			if (rotGroup.getRotation(i).getDirection() == 1) {
				axisAngle = rotGroup.getRotation(i).getAxisAngle();
				points[2] = new Point3d(axisAngle.x, axisAngle.y, axisAngle.z);
				points[2].scale(1);
				break;
			}
		}
		
		if (points[2] == null) {
			Vector3d v =  new Vector3d(subunits.getCenters().get(0));
			v.normalize();
			v.scale(1);
			points[2] = new Point3d(v);
		}
		return points;
	}
	
	private static double getAxisDirection(RotationGroup rotGroup, Subunits subunits) {
		Point3d[] points = new Point3d[3];
		// add points on principal axis 10 A in each direction
		AxisAngle4d axisAngle = rotGroup.getRotation(0).getAxisAngle();
		points[0] = new Point3d(axisAngle.x, axisAngle.y, axisAngle.z);
		points[0].scale(20);
		points[1] = new Point3d(axisAngle.x, axisAngle.y, axisAngle.z);
		points[1].scale(-20);
		
		double dNterminal0 = 0;
		double dNterminal1 = 0;
		double dCterminal0 = 0;
		double dCterminal1 = 0;
		List<Point3d[]> traces = subunits.getTraces();
		for (Point3d[] trace: traces) {
			int midpoint = trace.length/2;

			for (int i = 0; i < midpoint; i++) {
				dNterminal0+= trace[i].distanceSquared(points[0]);
				dNterminal1+= trace[i].distanceSquared(points[1]);
			}
			for (int i = midpoint; i < trace.length; i++) {
				dCterminal0+= trace[i].distanceSquared(points[0]);
				dCterminal1+= trace[i].distanceSquared(points[1]);
			}
		}
		System.out.println("dNterminal1/dNterminal0: " + dNterminal1/dNterminal0);
		System.out.println("dCterminal1/dCterminal0: " + dCterminal1/dCterminal0);
		System.out.println("ratio: " + (dNterminal1+dCterminal1)/(dNterminal0+dNterminal0));
		double ratio = (dNterminal1+dCterminal1)/(dNterminal0+dNterminal0);
		double sign = Math.signum(ratio-1);
		System.out.println("sign: " + sign);
		return sign;
	}
	
	private FindQuarternarySymmetry getQuaternarySymmetry(Structure structure) {
//	    Chain Cluster seqGroup = new GlobalSequenceGrouper(structure, MIN_SEQUENCE_LENGTH);
//		List<Point3d[]> caCoords = seqGroup.getCalphaCoordinates();
//		List<Point3d[]> cbCoords = seqGroup.getCbetaCoordinates();
//		List<Integer> sequenceClusterIds = seqGroup.getSequenceClusterIds();
		FindQuarternarySymmetry finder = new FindQuarternarySymmetry(structure, parameters);
		return finder;
	}
	
	private List<Integer> calcMapping(List<Point3d> originalCoords, List<Point3d> transformedCoords) {
		double distanceThreshold = 100;
		List<Integer> permutation = new ArrayList<Integer>(transformedCoords.size());
		double sum = 0.0f;
		boolean[] flagged = new boolean[originalCoords.size()];
		Arrays.fill(flagged, false);

		for (Point3d t: transformedCoords) {
			int closest = -1;
			double minDist = Double.MAX_VALUE;

			for (int j = 0; j < originalCoords.size(); j++) {
				if (!flagged[j]) {
					double dist = t.distanceSquared(originalCoords.get(j));
					if (dist < minDist) {
						closest = j;
						minDist = dist;
					} 
				}
			}

			sum += minDist;
			
			if (closest == -1) {
				break;
			}
			
			flagged[closest] = true;
			System.out.println("Closest mapping: " + closest + " dist: " + minDist);
			permutation.add(closest);
		}
		double rmsd = Math.sqrt(sum / transformedCoords.size());
		

		if (rmsd > distanceThreshold || permutation.size() != transformedCoords.size()) {
			permutation.clear();
			return permutation;
		}

		// check uniqueness of indices
		Set<Integer> set = new HashSet<Integer>(permutation);

		// if size mismatch, clear permutation (its invalid)
		if (set.size() != originalCoords.size()) {
			//      	System.out.println("RotationSolver: getPermutation: duplicate members" + set.size());
			permutation.clear();
		}

		//	        System.out.println("RMSD: " + rmsd + " missed: " + Math.sqrt(missed) + " closest: " + missedI + " count: " + count);
		//       System.out.println("P1: " + permutation);
		//       System.out.println("P2: " + p2);
		
		System.out.println("Mapping RMSD: " + rmsd);
		return permutation;
	}
}
