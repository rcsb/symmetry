package org.biojava3.structure.align.symm.quarternary;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.FileConvert;

public class QuatSymmetryWriter {
	private Structure structure;

	public QuatSymmetryWriter(Structure structure) {
		this.structure = structure;
	}
	
	public void writeTransformedStructure(Matrix4d matrix, String fileName) {
		Structure s = getTransformedStructure(matrix);
		FileConvert f = new FileConvert(s);
		String pdbFile = f.toPDB();
		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileWriter(fileName));
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		out.print(pdbFile);
		out.close();
	}
	
	private Structure getTransformedStructure(Matrix4d matrix) {
		Structure s = structure.clone();
		int models = 1;
		if (structure.isBiologicalAssembly()) {
			models = s.nrModels();
		}

		for (int i = 0; i < models; i++) {	
			for (Chain c: s.getChains(i)) {
				for (Group g: c.getAtomGroups()) {
					for (Atom a: g.getAtoms()) {
						double[] coords = a.getCoords();
						Point3d p = new Point3d(coords);
						matrix.transform(p);
                        p.get(coords);
						a.setCoords(coords);
					}
				}
			}
		}
		return s;
	}
}
