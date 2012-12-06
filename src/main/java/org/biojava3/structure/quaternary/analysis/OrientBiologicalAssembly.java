package org.biojava3.structure.quaternary.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;


public class OrientBiologicalAssembly {
	private String fileName = "";
	private String outputDirectory = "";
	private boolean verbose = false;

	public OrientBiologicalAssembly (String fileName, String outputDirectory, boolean verbose) {		
		this.fileName = fileName;
		this.outputDirectory = outputDirectory;
		this.verbose = verbose;
	}

	public static void main(String[] args) {		
		System.out.println("OrientBiologicalAssembly V 0.5: Calculates 4x4 transformation matrix to align structure along highest symmetry axis");
		System.out.println();
		if (args.length < 2) {
			System.out.println("Usage: OrientBiologicalAssembly pdbFile outputDirectory [-verbose]");
			System.exit(-1);
		}
		if (!args[0].contains(".pdb")) {
			System.err.println("Input file must be a .pdb or a pdb.gz file");
			System.exit(-1);
		}
		boolean verbose = false;
		if (args.length == 3 && args[2].equals("-verbose")) {
			verbose = true;
		}
		
		
		OrientBiologicalAssembly orienter = new OrientBiologicalAssembly(args[0], args[1], verbose);
		orienter.run();
	}

	public void run() {
		Structure structure = readStructure();
		orient(structure);
	}

	private Structure readStructure() {
		// initialize the PDB_DIR env variable
		AtomCache cache = new AtomCache();
		
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		p.setLoadChemCompInfo(true);
		p.setAtomCaThreshold(Integer.MAX_VALUE);
		p.setAcceptedAtomNames(new String[]{" CA "});

		PDBFileReader pdbreader = new PDBFileReader();
		pdbreader.setFileParsingParameters(p);

		Structure structure = null;
		try { 
			structure = pdbreader.getStructure(fileName);
			structure.setBiologicalAssembly(isBioassembly());
		} catch (Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
		return structure;
	}
	
	private void orient(Structure structure) {	
		
		// initialize with default parameters
		QuatSymmetryParameters params = new QuatSymmetryParameters();
		params.setSequenceIdentityThreshold(0.3);
		params.setVerbose(verbose);
		
		System.out.println("Default parameters:");
		System.out.println(params);
		
		
		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry();
		calc.setBioAssembly(structure);
		calc.setParams(params);
				
		calc.orient();
		
		AxisTransformation at = calc.getAxisTransformation();
		
		String prefix = getBaseFileName();
		String bioassemblyId = getBioassemblyId();
		if (! bioassemblyId.isEmpty()) {
		    prefix += "_" + bioassemblyId;
		}
		
		System.out.println("Bioassembly id: " + getBioassemblyId());
		
		String outName = prefix + "_4x4transformation.txt";
		System.out.println("Writing 4x4 transformation to: " + outName);
		writeFile(outName, at.getTransformation().toString());
		
		outName = prefix + "_JmolAnimation.txt";
		System.out.println("Writing Jmol animation to: " + outName);
		writeFile(outName, calc.getScriptGenerator().playOrientations());

		System.out.println("Color by subunit         : " + calc.getScriptGenerator().colorBySubunit());
		System.out.println("Color by subunit length  : " + calc.getScriptGenerator().colorBySubunit().length());
		System.out.println("Color by sequence cluster: " + calc.getScriptGenerator().colorBySequenceCluster());
		System.out.println("Color by seq. clst. len  : " + calc.getScriptGenerator().colorBySequenceCluster().length());
		System.out.println("Color by symmetry        : " + calc.getScriptGenerator().colorBySymmetry());
		System.out.println("Color by symmetry length : " + calc.getScriptGenerator().colorBySymmetry().length());

		System.out.println("Draw axes                : " + calc.getScriptGenerator().drawAxes());
		System.out.println("Draw axes length         : " + calc.getScriptGenerator().drawAxes().length());
		System.out.println("Draw polyhedron          : " + calc.getScriptGenerator().drawPolyhedron());
		System.out.println("Draw polyhedron length   : " + calc.getScriptGenerator().drawPolyhedron().length());

		
		System.out.println("Default orientation      : " + calc.getScriptGenerator().setDefaultOrientation());
		System.out.println("Orientation count        : " + calc.getScriptGenerator().getOrientationCount());
		for (int i = 0; i <  calc.getScriptGenerator().getOrientationCount(); i++) {
			System.out.println("Orientation name " + i + "       : " + calc.getScriptGenerator().getOrientationName(i));
			System.out.println("Orientation " + i + "            : " + calc.getScriptGenerator().setOrientation(i));
		}
		
		// avoid memory leaks
		calc.destroy();
	}
	
	private String getBioassemblyId() {
		int first = fileName.indexOf(".pdb") + 4;
		int last = fileName.length();
		int index = fileName.indexOf(".gz");
		if (index > 0) {
			last = index;
		}
		return fileName.substring(first, last);
	}
	
	private boolean isBioassembly() {
		return !getBioassemblyId().isEmpty();
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
		if (lastChar.equals('/') || lastChar.equals('\\')) {
			return outputDirectory + name;
		} else {
			return outputDirectory + File.separatorChar + name;
		}
	}
}
