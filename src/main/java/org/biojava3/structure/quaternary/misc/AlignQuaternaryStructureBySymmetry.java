package org.biojava3.structure.quaternary.misc;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.quaternary.analysis.CalcBioAssemblySymmetry;
import org.biojava3.structure.quaternary.core.ChainClusterer;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.QuatSymmetryWriter;
import org.biojava3.structure.quaternary.core.SequenceAlignmentCluster;
import org.biojava3.structure.quaternary.geometry.SuperPosition;

public class AlignQuaternaryStructureBySymmetry {
	private double SEQUENCE_IDENTITY_THRESHOLD = 0.95;
	private static String PDB_PATH = "C:/Users/Peter/Documents/PDB/";
	private static String FILENAME = "C:/Users/Peter/Documents/QuatStructureComparison/20130125_100800_symm95_pairs.csv";
	private static QuatSymmetryParameters parameters = new QuatSymmetryParameters();
	private Map<SequenceAlignmentCluster,SequenceAlignmentCluster> clusterMap = new LinkedHashMap<SequenceAlignmentCluster,SequenceAlignmentCluster>();
    private ClusterAlignment calign = new ClusterAlignment();
	private List<String[]> csv = new ArrayList<String[]>();
	private String header = "pdbId1,pdbId2,baId1,baId2,stoichiometry,pointgroup,technique1,technique2,spacegroup1,spacegroup2,resolution1,resolution2,classification1,classfication2, rmsdTert, rmsdTot, rmsdQuat";	

    
    public static void main(String[] args) throws IOException {//

    	// Aspartate Transcarbamylase (ATCase)
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1Q95","1R0B","1R0C","1RAA","1RAB","1RAC"});
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1R0B","2FZC"}); // D3, seq. length mismatch
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2FZC","1Q95"}); // T vs. R-state, error
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2J0X","2J0W"}); // + C2/A2 T vs. R-state, rmsd 9.85
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2ZQY","2ZQZ"}); // + D2/A4 (a+)T vs. R-state, rmsd 4.16
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1RDY","1RDX"}); // + D2/A4 (a+) T vs. R-state, rmsd 5.86
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1FRP","1YYZ"}); // + D2/A4 T vs. R-state, mismatch, rmsd 4.11 (required lower seq. id setting)
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1M35","1W2M"}); // + D2/A4 0.40
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1A4S","1BPW"}); // + D2/A4 0.43
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1AT1","1AT1"}); // + D3/A6B6
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1AHU","1AHV"}); // + D4/A8 1.15
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1AUS","1RCO"}); // + D4/A8B8 0.91
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1DHN","2NM2"}); // + D4/A8 0.89
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1UPM","1UPP"}); // + D4/A8B8 1.80
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1HO1","1HO4"}); // + D4/A8 0.60
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1SOR","3M9I"}); // - D4/A8 ??
    	//		List<String> pdbIds = Arrays.asList(new String[]{"3RTC","3RTD"}); // + D4/A8 0.26
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2C14","2C16"}); // + D4/A8 0.30
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2FW1","2FWJ"}); // + D4/A8 0.15
    	//		List<String> pdbIds = Arrays.asList(new String[]{"3BE7","3DUG"}q); // + D4/A8 2.9
    	//		List<String> pdbIds = Arrays.asList(new String[]{"3A8E","3AJ1"}); // + D4/A8 1.66
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2Q0Q","2Q0S"}); // + D4/A8 0.34
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2OJW","2QC8"}); // + D4/A8 D5/A10 0.69
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1BSR","1R5C"}); // + C2/A2 1.13
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1ADL","1LIB"}); // + C2/A2 0.31
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2OIE","2OIG"}); // + C2/A2 0.76
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1A3W","1A3X"}); // + D2/A4 0.98
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2ONO","2ONP"}); // + D2/A4 0.38
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2P3N","2P3V"}); // + D2/A4 2.28
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1A0S","1A0T"}); // + C3/A3 0.23
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1DF4","1DF5"}); // + C3/A3 0.81
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1B77","1B8H"}); // + C3/A3 0.51
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2OKD","2OKE"}); // + C3/A3 0.37
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2J58","2W8I"}); // + C8/A8// 0.41
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2BL2","2CYD"}); // + C10/A10 0.16
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1CAU","1CAV"}); // + C3/A3B3 0.85
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1DGR","1DGW"}); // + C3/A3B3C3 0.48
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1BCF","1BFR"}); // + O/A24 0.36
    	//		List<String> pdbIds = Arrays.asList(new String[]{"2IHW","2II4"}); // + O/A24 0.55

    	//		List<String> pdbIds = Arrays.asList(new String[]{"1AQ3","1MVB"}); // + I/A180 0.14
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1A16","1N51"}); // (+) D2/A4 0.45 need to check all 3 C2 axis in two directions = 6 possible alignments
    	//		List<String> pdbIds = Arrays.asList(new String[]{"4HHB","3R5I"}); // + C2/A2B2 3.55
    	//		List<String> pdbIds = Arrays.asList(new String[]{"1R0B","4HHB"}); // invalid combination for testing

    	// Transthyretin
    	//	List<String> pdbIds = Arrays.asList(new String[]{"3KGT","1TTC"});

 
   	AlignQuaternaryStructureBySymmetry aligner = new AlignQuaternaryStructureBySymmetry();
  //  	double[] rmsds = aligner.align("1EVR", "1M5A");
//       	double[] rmsds = aligner.align("4HHB", "3R5I");
 //      	double[] rmsds = aligner.align("4HHB", "1G9V");
 //   	double[] rmsds = aligner.align("1NH0", "2QNP");
    	double[] rmsds = aligner.align("1EVR", 1,"1M5A", 1);
    	
  //  	AlignQuaternaryStructureBySymmetry align = new AlignQuaternaryStructureBySymmetry();
	//	align.run();
    }
    public void run() throws IOException {

    	String outfile = FILENAME;
		outfile = outfile.substring(0, outfile.length()-4);
		String errfile = outfile;
		outfile += "_align.csv";
		errfile += "_align.txt";
		System.out.println("Output file: " + outfile);
		System.out.println("Error file: " + errfile);
		
		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileWriter(outfile));
		} catch (IOException e) {
			System.exit(-1);
		}
		out.println(header);
		
    	readFile(FILENAME);
    	int count = 0;
    	for (String[] items: csv) {
    		if (count > 0) {
    			System.out.println(Arrays.toString(items));
    			String pdbId1 = items[0];
    			String pdbId2 = items[1];
    			int baId1 = Integer.valueOf(items[2]);
    			int baId2 = Integer.valueOf(items[3]);	
    			double[] rmsd = align(pdbId1, baId1, pdbId2, baId2);
    			System.out.println("RMSD: " + rmsd[0] + " - " + rmsd[1]);
    			for (String s: items) {
    				out.print(s);
    				out.print(",");
    			}
    			out.print(rmsd[0]);
    			out.print(",");
    			out.print(rmsd[1]);
    			out.print(",");
    			out.println(rmsd[1]-rmsd[0]);
    		}
    		count++;
    	}
    }
	
	public double[] align(String pdbId1, int baId1, String pdbId2, int baId2) {
		double[] bestRmsd = {Double.MAX_VALUE, Double.MAX_VALUE};
		Structure s1 = getStructure(pdbId1, baId1);
		if (s1 == null) {
			return bestRmsd;
		}

		Structure s2 = getStructure(pdbId2, baId2);
		if (s2 == null) {
			return bestRmsd;
		}
		
		return align(s1, s2);
	}
	
	public double[] align(Structure s1, Structure s2) {
		double[] bestRmsd = {Double.MAX_VALUE, Double.MAX_VALUE};
		
		CalcBioAssemblySymmetry calc1 = orient(s1);	
//	    List<List<Integer>> orbits1 = calc1.getAxisTransformation().getOrbits();	
//	
//		System.out.println("Subunits 1: ");
//		for (List<Integer> orbit: orbits1) {
//			System.out.println(orbit);
//		}
//		System.out.println("Orient1: " + calc1.getScriptGenerator().getDefaultOrientation());
//	
//		CalcBioAssemblySymmetry calc2 = orient(s2);
//		
//		List<List<Integer>> orbits2 = calc2.getAxisTransformation().getOrbits();	
//		System.out.println("Subunits 2: ");
//		for (List<Integer> orbit: orbits2) {
//			System.out.println(orbit);
//		}
//		System.out.println("Orient2: " + calc2.getScriptGenerator().getDefaultOrientation());
//		
		ChainClustererNew grouper1 = new ChainClustererNew(s1, parameters);
		grouper1.addStructure(s2);
		System.out.println("Clusters:         " + grouper1.toString());
		System.out.println("Cluster ids:      " + grouper1.getSequenceClusterIds());
		System.out.println("Chain ids         " + grouper1.getChainIdsInClusterOrder());
		System.out.println("Chain ids         " + grouper1.getChainIds());
		System.out.println("Model nos:        " + grouper1.getModelNumbersInClusterOrder());
		System.out.println("Model nos:        " + grouper1.getModelNumbers());
		System.out.println("Struct ids:       " + grouper1.getStructureIdsInClusterOrder());
		System.out.println("Struct ids:       " + grouper1.getStructureIds());
		
//		for (List<Integer> orbit: orbits1) {
//			System.out.println(orbit);
//			boolean check = checkOrbit(orbit, 0, grouper1.getSequenceClusterIds());
//			if (! check) {
//				System.out.println("Orbit belongs to multiple sequence clusters");
//				return bestRmsd;
//			}
//		};
//		
//
//		int offset = calc1.getSubunits().getSubunitCount();
//		for (List<Integer> orbit: orbits2) {
//			System.out.println(orbit);
//			boolean check = checkOrbit(orbit, offset, grouper1.getSequenceClusterIds());
//			if (! check) {
//				System.out.println("Orbit belongs to multiple sequence clusters");
//				return bestRmsd;
//			}
//		};

//		List<Integer> clusterIds = grouper1.getSequenceClusterIds();
//		// calculate subunit RMSD
//		double aveRmsd = 0;
//		int count = 0;
//		List<Point3d[]> cas = grouper1.getCalphaCoordinates();
//		List<Point3d> xAll = new ArrayList<Point3d>();
//		List<Point3d> yAll = new ArrayList<Point3d>();
//		for (int i = 0; i < orbits1.size(); i++) {
//			System.out.println("Orbit: " + i);
//			List<Integer> orbit1 = orbits1.get(i);
//			List<Integer> orbit2 = orbits2.get(i);
//			if (orbit1.size() != orbit2.size()) {
//				System.out.println("Size of orbits doesn't match");
//				return bestRmsd;
//			}
//			for (int j = 0; j < orbit1.size(); j++) {
//				System.out.println("Subunit: " + j);
//				int index1 = orbit1.get(j);
//				int index2 = orbit2.get(j) + offset;
//				if (clusterIds.get(index1) != clusterIds.get(index2)) {
//					System.out.println("Cluster index mismatch");
//					return bestRmsd;
//				}
//				System.out.println(("indices: " + index1 + " - " + index2));
//				Point3d[] x = SuperPosition.clonePoint3dArray(cas.get(index1));
//				Point3d[] y = SuperPosition.clonePoint3dArray(cas.get(index2));
//		        xAll.addAll(Arrays.asList(cas.get(index1)));
//		        System.out.println("x.len: " + x.length);
//		        System.out.println("y.len: " + y.length);
//		        yAll.addAll(Arrays.asList(cas.get(index2)));
//		        
//			    Matrix4d m = SuperPosition.superposeWithTranslation(x, y);
//			    double rmsd = SuperPosition.rmsd(x, y);
//			    aveRmsd += rmsd*rmsd*x.length;
//			    count += x.length;
//			    System.out.println("RMSD: " + rmsd);
//			}
//		}
//	    bestRmsd[0] = Math.sqrt(aveRmsd/count);
//	    
//	    // calculate overall RMSD
//		Point3d[] x = xAll.toArray(new Point3d[0]);
//		x = SuperPosition.clonePoint3dArray(x);
//		 System.out.println(Arrays.toString(x));
//		Point3d[] y = yAll.toArray(new Point3d[0]);
//		y = SuperPosition.clonePoint3dArray(y);
//        
//	    Matrix4d m = SuperPosition.superposeWithTranslation(x, y);
//	    double totalRmsd = SuperPosition.rmsd(x, y);
//	    bestRmsd[1] = totalRmsd;
//		System.out.println("RESULT: rmsd: " + Arrays.toString(bestRmsd));
		return bestRmsd;
	}
	
	private boolean checkOrbit(List<Integer> orbit, int offset, List<Integer> clusterIds) {
		int refClusterId = clusterIds.get(orbit.get(0)+offset);
		for (int subunit: orbit) {
			if (clusterIds.get(subunit+offset) != refClusterId) {
				return false;
			}
			
		}
		return true;
	}

	private CalcBioAssemblySymmetry orient(Structure structure) {
		CalcBioAssemblySymmetry calc1 = new CalcBioAssemblySymmetry(structure, new QuatSymmetryParameters());
		calc1.orient();
//		Matrix4d transformation = calc1.getAxisTransformation().getTransformation();
//		transformStructure (structure, transformation);
		return calc1;
	}
	
	private void transformStructure(Structure structure, Matrix4d transformation) {
//		Structure s = structure.clone();
		int models = 1;
		if (structure.isBiologicalAssembly()) {
			models = structure.nrModels();
		}

		for (int i = 0; i < models; i++) {	
			for (Chain c: structure.getChains(i)) {
				for (Group g: c.getAtomGroups()) {
					for (Atom a: g.getAtoms()) {
						double[] coords = a.getCoords();
						Point3d p = new Point3d(coords);
						transformation.transform(p);
                        p.get(coords);
						a.setCoords(coords);
					}
				}
			}
		}
	}
	
	private Structure getStructure(String pdbId, int baId) {
		AtomCache cache = new AtomCache();
		cache.setAutoFetch(true);
		cache.setPath(PDB_PATH);	
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		//p.setParseCAOnly(true);
		//p.setMaxAtoms(50000000);
		
		p.setAtomCaThreshold(Integer.MAX_VALUE);
	//	p.setAcceptedAtomNames(new String[]{" CA ", " CB "});
		cache.setFileParsingParams(p);
			
		PrintWriter error = null;
		try {
			error = new PrintWriter(new FileWriter(PDB_PATH + "error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		Structure structure = null;
		try {
			structure = cache.getBiologicalAssembly(pdbId, baId, true);
			//				structure = cache. getStructure(pdbId);

		} catch (StructureException e) {
			e.printStackTrace();
			error.println(pdbId + "------------------------------------");
			error.println(e.getMessage());
			error.flush();
		} catch (IOException e) {
			e.printStackTrace();
			error.println(pdbId + "------------------------------------");
			error.println(e.getMessage());
			error.flush();
		}
		
		error.close();
		return structure;
	}
	
	private void readFile(String fileName) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] tokens = line.split(",");
			csv.add(tokens);
		}
		reader.close();
	}
}
