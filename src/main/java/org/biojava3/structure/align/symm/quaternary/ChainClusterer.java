package org.biojava3.structure.align.symm.quaternary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;

public class ChainClusterer  {
	private Structure structure = null;
	private QuatSymmetryParameters parameters = null;
	
	private List<Atom[]> caUnaligned = new ArrayList<Atom[]>();
	private List<String> chainIds = new ArrayList<String>();
	private List<String> sequences = new ArrayList<String>();
	private List<Atom[]> caAligned = new ArrayList<Atom[]>();
	private List<Atom[]> cbAligned = new ArrayList<Atom[]>();
	private List<Point3d[]> caCoords = new ArrayList<Point3d[]>();
	private List<Point3d[]> cbCoords = new ArrayList<Point3d[]>();

	List<SequenceAlignmentCluster> seqClusters = new ArrayList<SequenceAlignmentCluster>();
	
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
	
	public List<Point3d[]> getCbetaCoordinates() {
        run();
		return cbCoords;
	}
	
	public List<Atom[]> getCalphaTraces() {
		run();
		return caAligned;
	}
	
	public List<Atom[]> getCbetaTraces() {
		run();
		return cbAligned;
	}
	
	public boolean isHomomeric() {
		run();
		return seqClusters.size() == 1;
	}
	
	public int getMultiplicity() {
		run();
		return seqClusters.get(seqClusters.size()-1).getSequenceCount();
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
	
	public String getCompositionFormula() {
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
	
	private void run() {
		if (modified) {
			extractProteinChains();
			calcSequenceClusters();
			calcAlignedSequences();
			createCalphaTraces();
			createCbetaTraces();
			modified = false;
		}
	}
	
	private void extractProteinChains() {
		ProteinChainExtractor extractor = new ProteinChainExtractor(structure,  parameters.getMinimumSequenceLength());
		caUnaligned = extractor.getCalphaTraces();
		chainIds  = extractor.getChainIds();
		sequences = extractor.getSequences();
	    System.out.println("ChainClusterer: " + caUnaligned.size() + " " + chainIds);
        System.out.println("C alphas: ");
        for (Atom[] atoms: caUnaligned) {
        	System.out.println("Atoms: " + atoms.length);
        }
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
            UniqueSequenceList seqList = new UniqueSequenceList(caUnaligned.get(i), chainIds.get(i), sequences.get(i));
            SequenceAlignmentCluster seqCluster = new SequenceAlignmentCluster(parameters);
            seqCluster.addUniqueSequenceList(seqList);	
            seqClusters.add(seqCluster);
			
            for (int j = i + 1; j < caUnaligned.size(); j++) {
            	if (processed[j]) {
            		continue;
            	}
            	for (SequenceAlignmentCluster c: seqClusters) {
            		// add to existing sequence cluster if there is a match
            		if (c.isSequenceMatch(sequences.get(j))) {
            			if (c.addChain(caUnaligned.get(j), chainIds.get(j), sequences.get(j))) {
            				processed[j] = true;
            				break;
            			}
            		}
            	} 
            }

		}
		sortSequenceClustersBySize(seqClusters);
	}
	
	private void calcAlignedSequences() {
		caAligned = new ArrayList<Atom[]>();
		cbAligned = new ArrayList<Atom[]>();
		for (SequenceAlignmentCluster cluster: seqClusters) {
			caAligned.addAll(cluster.getAlignedCalphaAtoms());	
			cbAligned.addAll(cluster.getAlignedCBetaAtoms());
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
	
	private void createCbetaTraces() {
		for (Atom[] atoms: cbAligned) {
			Point3d[] trace = new Point3d[atoms.length];
			for (int j = 0; j < atoms.length; j++) {
				trace[j] = new Point3d(atoms[j].getCoords());
			}
			cbCoords.add(trace);
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
