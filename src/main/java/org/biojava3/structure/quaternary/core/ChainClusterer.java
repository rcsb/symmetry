package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;

public class ChainClusterer  {
	private Structure structure = null;
	private QuatSymmetryParameters parameters = null;
	
	private List<Atom[]> caUnaligned = new ArrayList<Atom[]>();
	private List<String> chainIds = new ArrayList<String>();
	private List<Integer> modelNumbers = new ArrayList<Integer>();
	private List<String> sequences = new ArrayList<String>();
	private List<Atom[]> caAligned = new ArrayList<Atom[]>();
	private List<Point3d[]> caCoords = new ArrayList<Point3d[]>();

	private List<SequenceAlignmentCluster> seqClusters = new ArrayList<SequenceAlignmentCluster>();
	
	private boolean modified = true;

	public ChainClusterer(Structure structure, QuatSymmetryParameters parameters) {
		this.structure = structure;
		this.parameters = parameters;
		modified = true;
	}	
	
	public List<Point3d[]> getCalphaCoordinates() {
        run();
		return caCoords;
	}
	
	public List<Atom[]> getCalphaTraces() {
		run();
		return caAligned;
	}
	
	public List<String> getChainIds() {
		run();
		List<String> chainIdList = new ArrayList<String>();

		for (int i = 0; i < seqClusters.size(); i++) {
	        SequenceAlignmentCluster cluster = seqClusters.get(i);
	        for (String chainId: cluster.getChainIds()) {
	        	chainIdList.add(chainId);
	        }
		}
		return chainIdList;
	}
	
	
	public List<Integer> getModelNumbers() {
		run();
		List<Integer> modNumbers = new ArrayList<Integer>();

		for (int i = 0; i < seqClusters.size(); i++) {
	        SequenceAlignmentCluster cluster = seqClusters.get(i);
	        for (Integer number: cluster.getModelNumbers()) {
	        	modNumbers.add(number);
	        }
		}
		return modNumbers;
	}
	
	public String getStoichiometry() {
		run();
		StringBuilder formula = new StringBuilder();
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

		for (int i = 0; i < seqClusters.size(); i++) {
			String c = "?";
			if (i < alpha.length()) {
				c = alpha.substring(i, i+1);
			}
			formula.append(c);
			int multiplier = seqClusters.get(i).getSequenceCount();
			if (multiplier > 1) {
				formula.append(multiplier);
			}
		}
		return formula.toString();
	}

	public List<Integer> getFolds() {
		run();
		List<Integer> denominators = new ArrayList<Integer>();
        Set<Integer> nominators = new TreeSet<Integer>();
		int nChains = caCoords.size();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			nominators.add(seqCount);
		}
		
		// find common denominators
		for (int d = 1; d <= nChains; d++) {

			int count = 0;
			for (Iterator<Integer> iter = nominators.iterator(); iter.hasNext();) {
				if (iter.next() % d == 0) {
					count++;
				}
			}
			if (count == nominators.size()) {
				denominators.add(d);
			}
		}
		
		return denominators;
	}
	
	public List<Integer> getSequenceClusterIds() {
		run();
		List<Integer> list = new ArrayList<Integer>();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			for (int i = 0; i < seqCount; i++) {
				list.add(id);
			}
		}
		return list;
	}
	
	
	public int getSequenceClusterCount() {
		run();
		return seqClusters.size();
	}
	
	public List<Boolean> getPseudoStoichiometry() {
		run();
		List<Boolean> list = new ArrayList<Boolean>();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			Boolean pseudo = seqClusters.get(id).isPseudoStoichiometric();
			for (int i = 0; i < seqCount; i++) {
				list.add(pseudo);
			}
		}
		return list;
	}
	
	public List<Double> getMinSequenceIdentity() {
		run();
		List<Double> list = new ArrayList<Double>();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			double minSequenceIdentity = seqClusters.get(id).getMinSequenceIdentity();
			for (int i = 0; i < seqCount; i++) {
				list.add(minSequenceIdentity);
			}
		}
		return list;
	}
	
	public List<Double> getMaxSequenceIdentity() {
		run();
		List<Double> list = new ArrayList<Double>();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			double maxSequenceIdentity = seqClusters.get(id).getMaxSequenceIdentity();
			for (int i = 0; i < seqCount; i++) {
				list.add(maxSequenceIdentity);
			}
		}
		return list;
	}
	
	private void run() {
		if (modified) {
			modified = false;
			
			extractProteinChains();
			if (caUnaligned.size() == 0) {
				return;
			}		
			
			calcSequenceClusters();
			mergeSequenceClusters();		
			calcAlignedSequences();
			createCalphaTraces();
		}
	}
	
	private void extractProteinChains() {
		ProteinChainExtractor extractor = new ProteinChainExtractor(structure,  parameters);
		caUnaligned = extractor.getCalphaTraces();
		chainIds  = extractor.getChainIds();
		sequences = extractor.getSequences();
		modelNumbers = extractor.getModelNumbers();
	}
	
	private void calcSequenceClusters() {
		boolean[] processed = new boolean[caUnaligned.size()];
		Arrays.fill(processed, false);
	
		for (int i = 0; i < caUnaligned.size(); i++) {
			if (processed[i]) {
				continue;
			}
			processed[i] = true;
			// create new sequence cluster
            UniqueSequenceList seqList = new UniqueSequenceList(caUnaligned.get(i), chainIds.get(i), modelNumbers.get(i), 0, sequences.get(i));
            SequenceAlignmentCluster seqCluster = new SequenceAlignmentCluster(parameters);
            seqCluster.addUniqueSequenceList(seqList);	
            seqClusters.add(seqCluster);
			
            for (int j = i + 1; j < caUnaligned.size(); j++) {
            	if (processed[j]) {
            		continue;
            	}
            	for (SequenceAlignmentCluster c: seqClusters) {
            			if (c.identityMatch(caUnaligned.get(j), chainIds.get(j), modelNumbers.get(j), 0, sequences.get(j))) {
            				processed[j] = true;
            				System.out.println("found identity match: " + i + " - " + j);
            				break;
            			}
            	} 
            }
		}
		sortSequenceClustersBySize(seqClusters);
	}
	
	private void mergeSequenceClusters() {
		boolean[] merged = new boolean[seqClusters.size()];
		Arrays.fill(merged, false);

		for (int i = 0, n = seqClusters.size(); i < n-1; i++) {
			if (! merged[i]) {
				SequenceAlignmentCluster c1 = seqClusters.get(i);
				for (int j = i + 1; j < n; j++) {
					SequenceAlignmentCluster c2 = seqClusters.get(j);
					int[][][] alignment = c1.alignClustersByStructure(c2);
					if (alignment != null) {
						merged[j] = true;
						System.out.println("Merged cluster: " + j + " -> " + i);
						mergeCluster(c1, c2, alignment);
					}
					//				System.out.println("Cluster strutural overlap: " + i + " - " + j + ": " + overlap);
				}
			}
		}
		for (int i = seqClusters.size()-1; i > 0; i--) {
			if (merged[i]) {
				System.out.println("removing merged cluster: " + i);
				seqClusters.remove(i);
			}
		}
		sortSequenceClustersBySize(seqClusters);
	}
	
	private void mergeCluster(SequenceAlignmentCluster c1,
			SequenceAlignmentCluster c2, int[][][] alignment) {

		for (UniqueSequenceList u: c2.getUniqueSequenceList()) {
			// set new sequence alignment
			List<Integer> align1 = new ArrayList<Integer>(alignment[0][0].length);
			for (Integer a1: alignment[0][0]) {
				align1.add(a1);
			}
			u.setAlignment1(align1);

			List<Integer> align2 = new ArrayList<Integer>(alignment[0][1].length);
			for (Integer a2: alignment[0][1]) {
				align2.add(a2);
			}
			u.setAlignment2(align2);	
			c1.addUniqueSequenceList(u);
		}
	}

	private void calcAlignedSequences() {
		caAligned = new ArrayList<Atom[]>();
		for (SequenceAlignmentCluster cluster: seqClusters) {
			caAligned.addAll(cluster.getAlignedCalphaAtoms());	
		}
	}
	
	private void createCalphaTraces() {
		for (Atom[] atoms: caAligned) {
			Point3d[] trace = new Point3d[atoms.length];
			for (int j = 0; j < atoms.length; j++) {
				trace[j] = new Point3d(atoms[j].getCoords());
			}
			caCoords.add(trace);
		}
	}

	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("Sequence alignment clusters: " + seqClusters.size());
		builder.append("\n");
		for (SequenceAlignmentCluster s: seqClusters) {
			builder.append("# seq: ");
			builder.append(s.getSequenceCount());
			builder.append(" alignment length: ");
			builder.append(s.getSequenceAlignmentLength());
			builder.append("\n");
		}
		return builder.toString();
	}
	
	public void sortSequenceClustersBySize(List<SequenceAlignmentCluster> clusters) {
		Collections.sort(clusters, new Comparator<SequenceAlignmentCluster>() {
			public int compare(SequenceAlignmentCluster c1, SequenceAlignmentCluster c2) {
				int sign = Math.round(Math.signum(c2.getSequenceCount() - c1.getSequenceCount()));
				if (sign != 0) {
					return sign;
				}
				return Math.round(Math.signum(c2.getSequenceAlignmentLength() - c1.getSequenceAlignmentLength()));
			}
		});
	}
}
