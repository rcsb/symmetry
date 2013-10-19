package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import javax.vecmath.Point3d;

public class SymmetryDeviation {
	private Subunits subunits = null;
	private List<List<Integer>> permutations = new ArrayList<List<Integer>>();
	private List<Integer> folds = new ArrayList<Integer>();

	public SymmetryDeviation(Subunits subunits, RotationGroup rotationGroup) {
		this.subunits = subunits;
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
			permutations.add(rotationGroup.getRotation(i).getPermutation());
			folds.add(rotationGroup.getRotation(i).getFold());
		}
	}
	
	public SymmetryDeviation(Subunits subunits, HelixLayers helixLayers) {
		this.subunits = subunits;
		permutations = new ArrayList<List<Integer>>();
		for (int i = 0; i < helixLayers.size(); i++) {
			permutations.add(helixLayers.getHelix(i).getPermutation());
			folds.add(1);
		}
	}
	
	public double getSymmetryDeviation() {
		List<List<Integer>> triplets = new ArrayList<List<Integer>>();
		List<Point3d[]> traces = subunits.getTraces();
		double sumTotal = 0;
		int countTotal = 0;
		double sumIsologous = 0;
		int countIsologous = 0;
		
		for (int i = 0; i < permutations.size(); i++) {
			List<Integer> permutation = permutations.get(i);
			int fold = folds.get(i);
			
			List<List<Integer>> layerLines = calcLayerLines(permutation);
			for (List<Integer> cycle: layerLines) {
				if (cycle.size() < 3) {
					continue;
				}
				System.out.println("Cycle: " + cycle);
				if (cycle.size() > 3 && cycle.get(0) == cycle.get(cycle.size()-1)) {
					cycle.remove(cycle.size()-1);
				}
				for (int j = 0; j < cycle.size()-2; j++) {
					List<Integer> triplet = new ArrayList<Integer>(3);
					triplet.add(cycle.get(j));
					triplet.add(cycle.get(j+1));
					triplet.add(cycle.get(j+2));
					if (triplet.get(0) > triplet.get(2)) {
						Collections.reverse(triplet);
					}
					if (triplets.contains(triplet)) {
						continue;
					}
					triplets.add(triplet);
					System.out.println("Comparing: " + triplet.get(0) + "-" + triplet.get(1) + " - " + triplet.get(1) + "-" + triplet.get(2));
					Point3d[] t0 = traces.get(triplet.get(0));
					Point3d[] t1 = traces.get(triplet.get(1));
					Point3d[] t2 = traces.get(triplet.get(2));

					for (int k = 0; k < t0.length; k++) {
						for (int l = 0; l < t0.length; l++) {
							double d1 = t0[k].distance(t1[l]);
							double d2 = t1[k].distance(t2[l]);			
							sumTotal += Math.abs(d1-d2);
							countTotal++;
							if (triplet.get(0) == triplet.get(2)) {
								sumIsologous += Math.abs(d1-d2);
								countIsologous ++;
							}
						}
					}
				}		
			}
		}
		double sumHeterologous = sumTotal - sumIsologous;
		int countHeterologous = countTotal - countIsologous;
		
		System.out.println("SymmetryDeviation total    : " + sumTotal/countTotal);
		System.out.println("SymmetryDeviation isologous: " + ((countIsologous == 0) ? 0 : sumIsologous/countIsologous));
		System.out.println("SymmetryDeviation heterologous: " + ((countHeterologous == 0) ? 0 : sumHeterologous/(countHeterologous)));
	
		return sumTotal/countTotal;
	}
	
	private List<List<Integer>> calcLayerLines(List<Integer> permutation) {
		List<List<Integer>> layerLines = new ArrayList<List<Integer>>();		
		createLineSegments(permutation, layerLines);		
		int count = layerLines.size();
		
		// iteratively join line segments
		do {
			count = layerLines.size();
			joinLineSegments(layerLines);
			// after joining line segments, get rid of the empty line segments left behind
			trimEmptyLineSegments(layerLines);
		} while (layerLines.size() < count);
		
		return layerLines;
	}

	private void createLineSegments(List<Integer> permutation,
			List<List<Integer>> layerLines) {
		for (int i = 0; i < permutation.size(); i++) {
			if (permutation.get(i) != -1 ) {
				List<Integer> lineSegment = new ArrayList<Integer>();
				lineSegment.add(i);
				lineSegment.add(permutation.get(i));
				layerLines.add(lineSegment);
			}
		}
	}
	
	private void joinLineSegments(List<List<Integer>> layerLines) {
		for (int i = 0; i < layerLines.size()-1; i++) {
			List<Integer> lineSegmentI = layerLines.get(i);
			if (! lineSegmentI.isEmpty()) {
				for (int j = i + 1; j < layerLines.size(); j++) {
					List<Integer> lineSegmentJ = layerLines.get(j);
					if (! lineSegmentJ.isEmpty()) {
						if (lineSegmentI.get(lineSegmentI.size()-1).equals(lineSegmentJ.get(0))) {
//							System.out.println("join right: " + lineSegmentI + " - " + lineSegmentJ);
							lineSegmentI.addAll(lineSegmentJ.subList(1,  lineSegmentJ.size()));
//							System.out.println("joned segment: " + lineSegmentI);
							lineSegmentJ.clear();		
						} else if ((lineSegmentI.get(0).equals(lineSegmentJ.get(lineSegmentJ.size()-1)))) {
							lineSegmentI.addAll(0, lineSegmentJ.subList(0,  lineSegmentJ.size()-1));
//							System.out.println("join left: " + lineSegmentJ + " - " + lineSegmentI);
//							System.out.println("joned segment: " + lineSegmentI);
							lineSegmentJ.clear();
						}
					}
				}
			}
		}
	}
	
	private void trimEmptyLineSegments(List<List<Integer>> layerLines) {
		for (Iterator<List<Integer>> iter = layerLines.iterator(); iter.hasNext();) {
			if (iter.next().isEmpty()) {
				iter.remove();
			}
		}
	}
}
