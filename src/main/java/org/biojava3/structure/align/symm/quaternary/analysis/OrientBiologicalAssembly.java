package org.biojava3.structure.align.symm.quaternary.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.FindQuarternarySymmetry;
import org.biojava3.structure.align.symm.quaternary.QuatSymmetryParameters;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;
import org.biojava3.structure.align.symm.quaternary.Subunits;

public class OrientBiologicalAssembly {
	private static final int MIN_SEQUENCE_LENGTH = 24;
	private static final double SEQUENCE_IDENTITY_THRESHOLD = 0.30;
	private static final double ALIGNMENT_FRACTION_THRESHOLD = 0.9;
	private static final double RMSD_THRESHOLD = 5.0;

	public OrientBiologicalAssembly () {
	}

	public static void main(String[] args) {
		Structure structure = readStructure(args[0]);
		Matrix4d matrix = orient(structure);
		outputRotation(args[1], matrix);
	}
	
	private static Structure readStructure(String filename) {
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		p.setLoadChemCompInfo(true);
		p.setAtomCaThreshold(Integer.MAX_VALUE);
		p.setAcceptedAtomNames(new String[]{" CA "});
//		p.setAcceptedAtomNames(new String[]{" CA ", " CB "});

		PDBFileReader pdbreader = new PDBFileReader();
		pdbreader.setFileParsingParameters(p);
		Structure structure = null;
		try{ 
			structure = pdbreader.getStructure(filename);
			structure.setBiologicalAssembly(true);

		} catch (Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
		return structure;
	}

	private static Matrix4d orient(Structure structure) {	
		QuatSymmetryParameters params = new QuatSymmetryParameters();
		params.setMinimumSequenceLength(MIN_SEQUENCE_LENGTH);
		params.setSequenceIdentityThreshold(SEQUENCE_IDENTITY_THRESHOLD);
		params.setAlignmentFractionThreshold(ALIGNMENT_FRACTION_THRESHOLD);
		params.setRmsdThreshold(RMSD_THRESHOLD);

		FindQuarternarySymmetry finder = new FindQuarternarySymmetry(structure, params);
		System.out.println("Formula: " + finder.getCompositionFormula());

		// return identity matrix if no protein chains are found
		if (finder.getChainCount() == 0) {
			Matrix4d matrix = new Matrix4d();
			matrix.setIdentity();
			return matrix;
		}

		RotationGroup rotationGroup = finder.getRotationGroup();	
		System.out.println("Point group: " + rotationGroup.getPointGroup());		
		System.out.println("Symmetry RMSD: " + (float) rotationGroup.getAverageTraceRmsd());

		Subunits subunits = finder.getSubunits();
		List<String> chainIds = finder.getChainIds();
		AxisTransformation at = new AxisTransformation(subunits, rotationGroup, chainIds);
		Matrix4d matrix = at.getTransformation();
		return matrix;
	}

	private static void outputRotation(String filename, Matrix4d matrix) {
		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileWriter(filename));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		out.println("transformation");
		out.println(matrix);
		out.flush();
		out.close();
		System.out.println("writing transformation");
	}
}
