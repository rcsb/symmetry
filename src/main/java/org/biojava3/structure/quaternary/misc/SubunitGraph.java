package org.biojava3.structure.quaternary.misc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.vecmath.Matrix4d;
import javax.vecmath.Matrix4f;
import javax.vecmath.Point3d;
import javax.vecmath.Point3f;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava3.structure.quaternary.core.Subunits;

public class SubunitGraph {
	private Structure structure = null;
    private List<Chain> chains;
	private List<Point3d[]> traces = new ArrayList<Point3d[]>();
	private Subunits subunits = null;

//
//	public SubunitGraph(Structure structure) {
//		this.structure = structure;
//	}
//
//	private Graph<Integer> getProteinGraph() {
//		int n = cbTraces.size();
//
//		// add vertex for each chain center
//		Graph<Integer> graph = new SimpleGraph<Integer>(); 
//		for (int i = 0; i < n; i++) {
//			graph.addVertex(i);
//		}
//
//		// add edges if there are 10 or more contact of CB atoms
//		for (int i = 0; i < n - 1; i++) {
//			for (int j = i + 1; j < n; j++) {
//				if (calcContactNumber(cbTraces.get(i), cbTraces.get(j)) >= 10) {
//					graph.addEdge(i, j);
//				}
//			}
//		}
//
//		return graph;
//	}
//
//	private int calcContactNumber(List<Point3d> a, List<Point3d> b) {
//		int contacts = 0;
//		for (Point3d pa : a) {
//			for (Point3d pb : b) {
//				if (pa.distance(pb) < 12) {
//					contacts++;
//				}
//			}
//		}
//		return contacts;
//	}
//
//	private Point3d calcCentroid(List<Point3d> points) {
//		Atom[] dummies = new Atom[points.size()];
//		for (int i = 0; i < points.size(); i++) {
//			double[] c = new double[3];
//			points.get(i).get(c);
//			dummies[i] = new AtomImpl();
//			dummies[i].setCoords(c);
//		}
//		Atom cent = Calc.getCentroid(dummies);
//		Point3d centroid = new Point3d(cent.getCoords());
//		return centroid;
//	}
//
//	private void transformCoordinates(Matrix4d matrix, Point3d centerOfMass) {
//		Structure s = structure.clone();
//		int models = 1;
//		if (biologicalAssembly) {
//			models = s.nrModels();
//		}
//
//		for (int i = 0; i < models; i++) {	
//			for (Chain c: s.getChains(i)) {
//				for (Group g: c.getAtomGroups()) {
//					for (Atom a: g.getAtoms()) {
//						Point3d p = new Point3d(a.getCoords());
//						p.sub(centerOfMass);
//						matrix.transform(p, p);
//					}
//				}
//			}
//		}
//	}
}
