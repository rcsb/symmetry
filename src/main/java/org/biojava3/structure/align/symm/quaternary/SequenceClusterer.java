package org.biojava3.structure.align.symm.quaternary;

import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;

public interface SequenceClusterer {

	public abstract List<Point3d[]> getCalphaCoordinates();

	public abstract List<Point3d[]> getCbetaCoordinates();

	public abstract List<Atom[]> getCalphaTraces();

	public abstract List<Atom[]> getCbetaTraces();

	public abstract boolean isHomomeric();

	public abstract int getMultiplicity();

	public abstract List<String> getOrderedChainIDList();

	public abstract String getCompositionFormula();

	public abstract List<Integer> getSequenceClusterIds();

}