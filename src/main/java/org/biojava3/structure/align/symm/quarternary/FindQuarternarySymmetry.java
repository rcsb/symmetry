package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.List;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;

public class FindQuarternarySymmetry {
	private List<Point3d[]> caCoords = null;
	private List<Point3d[]> cbCoords = null;
	private List<Integer> sequenceClusterIds = null;
	private boolean pseudoSymmetryAllowed = false;
	private Subunits subunits = null;
	private RotationGroup symmetryOperations = null;
	private String method = "";


	public FindQuarternarySymmetry(List<Point3d[]> caCoords, List<Point3d[]> cbCoords, List<Integer> sequenceClusterIds) {
		this.caCoords = caCoords;
		this.cbCoords = cbCoords;
		this.sequenceClusterIds = sequenceClusterIds;
	}
	
	public void setPseudoSymmetryAllowed(boolean pseudoSymmetryAllowed) {
		this.pseudoSymmetryAllowed = pseudoSymmetryAllowed;
	}

	public RotationGroup getRotationGroup() {
        long t1 = System.nanoTime();
        subunits = new Subunits(caCoords, cbCoords, sequenceClusterIds);
 //       System.out.println("Chains: " + traces.size());
        QuatSymmetryPerceptor perceptor = new QuatSymmetryPerceptor(subunits);
        perceptor.setPseudoSymmetryAllowed(pseudoSymmetryAllowed);
        symmetryOperations = perceptor.getSymmetryOperations();
        method = perceptor.getMethod();

 //       System.out.println("--- SymmetryOperations ---");;
 //       System.out.println(symmetryOperations);
//        System.out.println(symmetryOperations.getPointGroup());
//        createOutput();
        long t2 = System.nanoTime();
 //       System.out.println("Total Time: " + Math.round((t2 - t1) / 1000000) + " msec.");
        return symmetryOperations;
	}
	
	public Subunits getSubunits() {
		return subunits;
	}
	
	public int getChainCount() {
		return caCoords.size();
	}
	
	public String getMethod() {
		return method;
	}
}
