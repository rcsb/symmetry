package org.biojava3.structure.align.symm.quaternary;

import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;

public interface SequenceClusterer {

	public abstract void setStructure(Structure structure);

	/**
	 * @param minSequenceLength the minSequenceLength to set
	 */
	public abstract void setMinSequenceLength(int minSequenceLength);

	public abstract List<Point3d[]> getCalphaCoordinates();

	public abstract List<Point3d[]> getCbetaCoordinates();

	public abstract List<Atom[]> getCalphaTraces();

	public abstract List<Atom[]> getCbetaTraces();

	public abstract List<Chain> getChains();

	public abstract boolean isSequenceNumberedCorrectly();

	public abstract boolean isUnknownSequence();

	public abstract List<String[]> getSequences();

	public abstract boolean isHomomeric();

	public abstract int getMultiplicity();

	public abstract List<String> getOrderedChainIDList();

	public abstract String getCompositionFormula();

	public abstract List<List<Integer>> getSequenceCluster100();

	public abstract List<Integer> getSequenceClusterIds();

	public abstract List<Integer> getSequenceIds();

	public abstract void sortClusterBySize(List<List<Integer>> clusters);

}