package org.biojava3.structure.align.quaternary;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;

public class AlignTest {
	private static String PDB_PATH = "C:/PDB/";
	private Map<String,List<String>> structureGroups = new HashMap<String,List<String>>();
	private Map<String, String[]> csv = new HashMap<String,String[]>();

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
		
		AligQuaternaryStructure aligner = new AligQuaternaryStructure();
		double[] rmsds = aligner.align("1EVR", "1M5A");
		
	//	AlignTest test = new AlignTest();
	//	test.readFile(fileName);
	//	test.createStructureGroups();
	//	test.run();

	}
	
	private void run() {
		PrintWriter out = null;

		try {
			out = new PrintWriter(new FileWriter(PDB_PATH + "quat_sym.csv"));
		} catch (IOException e) {
			System.exit(-1);
		}
		
		out.println("pdbId1,pdbId2,composition,signature,quatRmsd,tertRmsd,deltaRmsd,pointGroup,spaceGroup1,spaceGroup2,technique1,technique2,ligand1,ligands2,interactiontype1,interactiontype2");
  
		System.out.println("Number of clusters: " + structureGroups.size());
		
		int comparisons = 0;
		for (Entry<String,List<String>> entry: structureGroups.entrySet()) {
			List<String> pdbIds = entry.getValue();
			if (pdbIds.size() > 1) {
				comparisons += (pdbIds.size()* pdbIds.size() - pdbIds.size())/2;
			}
		}
		System.out.println("TOTAL comparisons: " + comparisons);
		int count = 1;
		for (Entry<String,List<String>> entry: structureGroups.entrySet()) {
			List<String> pdbIds = entry.getValue();
			if (pdbIds.size() > 1) {
				System.out.println("+++++++++++ Aligning: " + count + "/" + comparisons + ": " + entry.getKey() + "+++++++++++++++++");
				for (int i = 0, n = pdbIds.size(); i < n-1; i++) {
					String pdbIdI = pdbIds.get(i);
					System.out.println("pdbIdI: " + pdbIdI);
					Structure s1 = getStructure(pdbIdI);
					if (s1 == null) {
						continue;
					}
			//		PDBCrystallographicInfo info1 = s1.getCrystallographicInfo();
					PDBHeader header1 = s1.getPDBHeader();
					String t1 = header1.getTechnique();
				//	float r1 = header1.getResolution(); resolution is missing in BU files
//					System.out.println("Classification: " + header1.getClassification());
//					System.out.println("Description: " + header1.getDescription());
//					System.out.println("Title: " + header1.getTitle());
			
//					if (! pdbIdI.equals("1HB9")) continue;
					for (int j = i+1; j < n; j++) {
						String pdbIdJ = pdbIds.get(j);
	//					if (! pdbIdJ.equals("1HB7")) continue;
						Structure s2 = getStructure(pdbIdJ);
						if (s2 == null) {
							continue;
						}
						System.out.println("pdbIdI: " + pdbIdJ);
			//			PDBCrystallographicInfo info2 = s2.getCrystallographicInfo();
						PDBHeader header2 = s2.getPDBHeader();
						String t2 = header2.getTechnique();
					//	float r2 = header2.getResolution();
						AligQuaternaryStructure aligner = new AligQuaternaryStructure();
				 //       double[] rmsds = aligner.align(pdbIdI, pdbIdJ);
				        double[] rmsds = aligner.align(s1, s2);
				        double deltaRmsd = rmsds[0]-rmsds[1];
				        out.println(pdbIdI + "," + pdbIdJ +"," + getComposition(pdbIdI) + "," + 
				        		entry.getKey() + "," + 
				        		rmsds[0] + "," + rmsds[1] + "," + deltaRmsd + "," +
				        		getPointGroup(pdbIdI) + "," +
	//			        		info1.getSpaceGroup() + "," + info2.getSpaceGroup() + "," +
				        		t1 + "," + t2 + "," +
				        		getLigands(pdbIdI) + "," +
				        		getLigands(pdbIdJ) + "," +
				                getInteractionType(pdbIdI) + "," +
		        	           	getInteractionType(pdbIdJ));
				        out.flush();
						count++;
					}
				}
			}

		}
		out.close();
	}
	
	private String getComposition(String pdbId) {
		return csv.get(pdbId)[1];
	}
	
	private String getPointGroup(String pdbId) {
		return csv.get(pdbId)[14];
	}
	private String getLigands(String pdbId) {
		String[] tokens = csv.get(pdbId);
        return tokens[28];
	}
	
	private String getInteractionType(String pdbId) {
		String[] tokens = csv.get(pdbId);
        return tokens[29];
	}
	private void readFile(String fileName) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] tokens = line.split(",");
			String pdbId = tokens[0];
			csv.put(pdbId, tokens);
		}
	}
	
	private void createStructureGroups() {
		for (Entry<String, String[]> entry: csv.entrySet()) {
			String pdbId = entry.getKey();
			String[] tokens = entry.getValue();
			String signature = tokens[2];
			String pointGroup = tokens[3];
			String key = signature + "_" + pointGroup;
		    if (structureGroups.containsKey(key)) {
		    	List<String> pdbIds = structureGroups.get(key);
		    	pdbIds.add(pdbId);
		    } else {
		    	List<String> list = new ArrayList<String>();
		    	list.add(pdbId);
		    	structureGroups.put(key, list);
		    }
		}
	}
	
	private Structure getStructure(String pdbId) {
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
			structure = cache.getBiologicalAssembly(pdbId, 1, true);
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
		return structure;
	}
}
