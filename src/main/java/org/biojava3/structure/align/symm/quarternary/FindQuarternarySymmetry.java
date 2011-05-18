package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;

public class FindQuarternarySymmetry {
	private Structure structure = null;
	private int minSequenceLength = 30;
	private int maxSequenceLengthDeviation = 20;

	/**
	 * @param minSequenceLength the minSequenceLength to set
	 */
	public void setMinSequenceLength(int minSequenceLength) {
		this.minSequenceLength = minSequenceLength;
	}

	/**
	 * @param maxSequenceLengthDeviation the maxSequenceLengthDeviation to set
	 */
	public void setMaxSequenceLengthDeviation(int maxSequenceLengthDeviation) {
		this.maxSequenceLengthDeviation = maxSequenceLengthDeviation;
	}
	
	public void run(Structure structure) {
		this.structure = structure;
		List<Atom> centroids = getChainCentroids();
		System.out.println(structure.getPDBCode() + " chains: " + centroids.size());
	}
	
	private List<Atom> getChainCentroids() {
		List<Atom> centroids = new ArrayList<Atom>();
		List<Atom[]> traces = getCAChainTraces();
		for (Atom[] trace: traces) {
			centroids.add(Calc.getCentroid(trace));
		}
		return centroids;
	}
	
	private List<Atom[]> getCAChainTraces() {
		System.out.println("getCAChainTraces");
		List<Atom[]> traces = new ArrayList<Atom[]>();

		for (int i = 0; i < structure.nrModels(); i++) {
			List<Chain> chains = structure.getChains(i);
		
			for (Chain c: chains) {
				Atom[] caAtoms = StructureTools.getAtomCAArray(c);
				if (caAtoms.length > minSequenceLength) {
					traces.add(StructureTools.getAtomCAArray(c));
				}
			}
		}
		
		return traces;
	}

}
