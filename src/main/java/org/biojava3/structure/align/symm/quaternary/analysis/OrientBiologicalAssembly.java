package org.biojava3.structure.align.symm.quaternary.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
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

	private String fileName = "";
	private String outputDirectory = "";

	public OrientBiologicalAssembly (String fileName, String outputDirectory) {
		this.fileName = fileName;
		this.outputDirectory = outputDirectory;
	}

	public static void main(String[] args) {
		System.out.println("Calculates 4x4 transformation matrix to align structure along highest symmetry axis");
		System.out.println();
		if (args.length != 2) {
			System.out.println("Usage: OrientBiologicalAssembly pdbFile outputDirectory");
			System.exit(-1);
		}
		OrientBiologicalAssembly orienter = new OrientBiologicalAssembly(args[0], args[1]);
		orienter.run();
	}

	public void run() {
		System.out.println("Default parameters:");
		System.out.println("Minimum protein sequence length: " + MIN_SEQUENCE_LENGTH);
		System.out.println("Sequence identity threshold    : " +  Math.round(SEQUENCE_IDENTITY_THRESHOLD*100) + "%");
		System.out.println("Alignment length threshold     : " +  Math.round(ALIGNMENT_FRACTION_THRESHOLD*100) + "%");
		System.out.println("Symmetry RMSD threshold        : " +  RMSD_THRESHOLD);
		System.out.println();

		Structure structure = readStructure(fileName);
		System.out.println("Protein chains used for alignment:");
		orient(structure);
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

	private void orient(Structure structure) {	
		QuatSymmetryParameters params = new QuatSymmetryParameters();
		params.setMinimumSequenceLength(MIN_SEQUENCE_LENGTH);
		params.setSequenceIdentityThreshold(SEQUENCE_IDENTITY_THRESHOLD);
		params.setAlignmentFractionThreshold(ALIGNMENT_FRACTION_THRESHOLD);
		params.setRmsdThreshold(RMSD_THRESHOLD);


		FindQuarternarySymmetry finder = new FindQuarternarySymmetry(structure, params);	


		RotationGroup rotationGroup = new RotationGroup();
		if (finder.getChainCount() > 0) {
			System.out.println();
			rotationGroup = finder.getRotationGroup();
			System.out.println("Results for " + Math.round(SEQUENCE_IDENTITY_THRESHOLD*100) + "% sequence identity threshold:");
			System.out.println("Stoichiometry: " + finder.getCompositionFormula());
			System.out.println("Point group  : " + rotationGroup.getPointGroup());		
			System.out.println("Symmetry RMSD: " + (float) rotationGroup.getAverageTraceRmsd());
		} 

		AxisTransformation at = new AxisTransformation(finder.getSubunits(), rotationGroup, finder.getChainIds());
	
		System.out.println();
		String outName = getBaseFileName() + "_4x4transformation.txt";
		System.out.println("Writing 4x4 transformation to: " + outName);
		writeFile(outName, at.getTransformation().toString());
		
		outName = getBaseFileName() + "_JmolTransformation.txt";
		System.out.println("Writing Jmol transformation to: " + outName);
		writeFile(outName, at.getJmolTransformation());
		
		outName = getBaseFileName() + "_JmolSymmetryAxes.txt";
		System.out.println("Writing Jmol symmetry axes to: " + outName);
		writeFile(outName, at.getJmolSymmetryAxes());
	}

	private static void writeFile(String fileName, String text) {
		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileWriter(fileName));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		out.println(text);
		out.flush();
		out.close();
	}

	private String getBaseFileName() {
		File f = new File(fileName);
		String name = f.getName();
		name = name.substring(0, name.indexOf('.'));
		Character lastChar = outputDirectory.charAt(outputDirectory.length()-1);
		if (lastChar.equals('/')) {
			return outputDirectory + name;
		} else {
			return outputDirectory + "\\" + name;
		}
	}
}
