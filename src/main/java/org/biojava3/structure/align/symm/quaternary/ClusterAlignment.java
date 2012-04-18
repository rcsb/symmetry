package org.biojava3.structure.align.symm.quaternary;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ClusterAlignment {
    private List<List<Integer>> clusters1 = new ArrayList<List<Integer>>();
    private List<List<Integer>> clusters2 = new ArrayList<List<Integer>>();
    private List<List<Integer>> alignments1 = new ArrayList<List<Integer>>();
    private List<List<Integer>> alignments2 = new ArrayList<List<Integer>>();
    private List<Double> sequenceIdentities = new ArrayList<Double>();
    
    public void add(List<Integer> cluster1, List<Integer> cluster2, List<Integer> alignment1, List<Integer> alignment2, double sequenceIdentity) {
    	clusters1.add(cluster1);
    	clusters2.add(cluster2);
//    	if (cluster1.size() != cluster2.size()) {
//    		System.out.println("ClusterAlignment: ERROR: cluster size mismatch: "+  cluster1.size() + " vs. " + cluster2.size());
//    	}
    	alignments1.add(alignment1);
    	alignments2.add(alignment2);
//    	if (alignment1.size() != alignment2.size()) {
//    		System.out.println("ClusterAlignment: ERROR: alignment size mismatch: "+  alignment1.size() + " vs. " + alignment2.size());
//    	}
    	sequenceIdentities.add(sequenceIdentity);
    }
    
    public int getAlignmentCount() {
    	return alignments1.size();
    }
    
    public List<Integer> getCluster1(int index) {
    	return clusters1.get(index);
    }
    
    public List<Integer> getCluster2(int index) {
    	return clusters2.get(index);
    }
    
    public List<Integer> getAlignment1(int index) {
    	return alignments1.get(index);
    }
    
    public List<Integer> getAlignment2(int index) {
    	return alignments2.get(index);
    }
//    public List<Integer> getAlignment1(int subunitId1) {
//    	for (List<Integer> cluster1: clusters1) {
//    		if (cluster1.contains(subunitId1)) {
//    			return alignments1.get(clusters1.indexOf(cluster1));
//    		}
//    	}
//    	return Collections.emptyList();
//    }
//
//    public List<Integer> getAlignment2(int subunitId2) {
//    	for (List<Integer> cluster2: clusters2) {
//    		if (cluster2.contains(subunitId2)) {
//    			return alignments2.get(clusters2.indexOf(cluster2));
//    		}
//    	}
//    	return Collections.emptyList();
//    }
}
