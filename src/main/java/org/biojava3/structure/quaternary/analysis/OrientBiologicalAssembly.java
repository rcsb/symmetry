package org.biojava3.structure.quaternary.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.vecmath.Matrix4d;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.ChemCompProvider;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.QuatSymmetryDetector;
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
		System.out.println("OrientBiologicalAssembly V " + CalcBioAssemblySymmetry.version + " - "  +  CalcBioAssemblySymmetry.build + " : Calculates 4x4 transformation matrix to align structure along highest symmetry axis");
		System.out.println();
		
		
//		AllChemCompProvider all = new AllChemCompProvider();
		ChemCompProvider all = new DownloadChemCompProvider();
		
		ChemCompGroupFactory.setChemCompProvider(all);
		
		// enforce the use of the version of vecmath that we bundle with this jar
		try {
			Class.forName("javax.vecmath.Point2i");
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		
		
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
		params.setVerbose(verbose);

		System.out.println("Default parameters:");
		System.out.println(params);

		long t1 = System.nanoTime();
		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry(structure, params);

		String prefix = getBaseFileName();
		String bioassemblyId = getBioassemblyId();
		if (! bioassemblyId.isEmpty()) {
			prefix += "_" + bioassemblyId;
		}

		
		QuatSymmetryDetector detector = calc.orient();
		
		boolean hasProtein = detector.hasProteinSubunits();
		
		long t2 = System.nanoTime();
		
		System.out.println("CalcBioAssemblySymmetry: " + (t2-t1)*0.000001 + " ms");
//
//		if (hasProtein) {		
//			
//			System.out.println("Bioassembly id: " + getBioassemblyId());
//			
//			System.out.println("Point group:             : " + calc.getRotationGroup().getPointGroup());
//			System.out.println("Subunit count            : " + calc.getSubunits().getSubunitCount());
//			System.out.println("Stoichiometry            : " + calc.getSubunits().getStoichiometry());
//			System.out.println("Color by subunit         : " + calc.getScriptGenerator().colorBySubunit());
//			System.out.println("Color by sequence cluster: " + calc.getScriptGenerator().colorBySequenceCluster());
//			System.out.println("Color by symmetry        : " + calc.getScriptGenerator().colorBySymmetry());
//
//			System.out.println("Draw axes                : " + calc.getScriptGenerator().drawAxes());
//			System.out.println("Draw polyhedron          : " + calc.getScriptGenerator().drawPolyhedron());
//
//			System.out.println("Zoom                     : " + calc.getScriptGenerator().getZoom());
//			System.out.println("Default orientation      : " + calc.getScriptGenerator().getDefaultOrientation());
//			
//			System.out.println("Orientation count        : " + calc.getScriptGenerator().getOrientationCount());
//			
//			for (int i = 0; i <  calc.getScriptGenerator().getOrientationCount(); i++) {
//				System.out.println("Orientation name " + i + "       : " + calc.getScriptGenerator().getOrientationName(i));
//				System.out.println("Orientation " + i + "            : " + calc.getScriptGenerator().getOrientation(i));
//				System.out.println("Orientation with zoom " + i + "  : " + calc.getScriptGenerator().getOrientationWithZoom(i));
//			}
//
//			String outName = prefix + "_4x4transformation.txt";
//			System.out.println("Writing 4x4 transformation to: " + outName);
//			AxisTransformation at = calc.getAxisTransformation();	
//			writeFile(outName, at.getTransformation().toString());
//
//			outName = prefix + "_JmolAnimation.txt";
//			System.out.println("Writing Jmol animation to: " + outName);
//			writeFile(outName, calc.getScriptGenerator().playOrientations());
//			
//		} else { 
//			System.out.println("Bioassembly id: " + getBioassemblyId());	
//			System.out.println("No protein chain found: returning identity matrix");
//			String outName = prefix + "_4x4transformation.txt";
//			System.out.println("Writing 4x4 transformation to: " + outName);
//			Matrix4d m = new Matrix4d();
//			m.setIdentity();
//			writeFile(outName, m.toString());
//		}
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
