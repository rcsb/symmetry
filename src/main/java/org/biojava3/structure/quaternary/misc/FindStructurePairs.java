package org.biojava3.structure.quaternary.misc;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava.bio.structure.PDBCrystallographicInfo;
import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.xtal.SpaceGroup;

public class FindStructurePairs {
	private static String PDB_PATH = "C:/Users/Peter/Documents/PDB/";
	private static String FILENAME = "C:/Users/Peter/Documents/QuatStructureComparison/20130125_100800_symm95.csv";
	// indices to data in input csv file
	private static int PDB_ID = 0;
	private static int BIOASSEMBLY = 1;
	private static int FORMULA = 2;
	private static int SIGNATURE95 = 5;
	private static int STOICHIOMETRY95 = 6;
	private static int POINT_GROUP = 9;
	
	private Map<String,List<String>> structureGroups = new HashMap<String,List<String>>();
	private Map<String, String[]> csv = new HashMap<String,String[]>();
	
	private String outfile = null;
	private String errfile = null;

	public static void main(String[] args) throws IOException {//
		FindStructurePairs test = new FindStructurePairs();
		test.readFile(FILENAME);
		test.createStructureGroups();
		test.run();

	}
	
	private void run() {	
		outfile = FILENAME;
		outfile = outfile.substring(0, outfile.length()-4);
		errfile = outfile;
		outfile += "_pairs.csv";
		errfile += "_pairs.txt";
		System.out.println("Output file: " + outfile);
		System.out.println("Error file: " + errfile);
		
		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileWriter(outfile));
		} catch (IOException e) {
			System.exit(-1);
		}
		
		PrintWriter error = null;
		try {
			error = new PrintWriter(new FileWriter(errfile));
		} catch (IOException e) {
			System.exit(-1);
		}
	    
		
		out.println("pdbId1,pdbId2,baId1,baId2,stoichiometry,pointgroup,technique1,technique2,spacegroup1,spacegroup2,resolution1,resolution2,classification1, classfication2");
  
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
			List<String> structureIds = entry.getValue();
			if (structureIds.size() > 1) {
				System.out.println("+++++++++++ Aligning: " + count + "/" + comparisons + ": " + entry.getKey() + "+++++++++++++++++");
				for (int i = 0, n = structureIds.size(); i < n-1; i++) {
					String pdbId1 = structureIds.get(i).substring(0, 4);
					String ba1 = structureIds.get(i).substring(4, 5);
					Structure s1 = getStructure(pdbId1, error);
					if (s1 == null) {
						continue;
					}
					PDBCrystallographicInfo info1 = s1.getCrystallographicInfo();
					SpaceGroup sp1 = info1.getSpaceGroup();
					PDBHeader header1 = s1.getPDBHeader();
					String t1 = header1.getTechnique();
					float r1 = header1.getResolution();
					String c1 = header1.getClassification();
		//			System.out.println("Classification: " + header1.getClassification());
		//			System.out.println("Description: " + header1.getDescription());
		//			System.out.println("Title: " + header1.getTitle());
			
					for (int j = i+1; j < n; j++) {
						String pdbId2 = structureIds.get(j).substring(0, 4);
						String ba2 = structureIds.get(j).substring(4, 5);
						Structure s2 = getStructure(pdbId2, error);
						if (s2 == null) {
							continue;
						}
						System.out.println("pdbIdJ: " + pdbId2);
						PDBCrystallographicInfo info2 = s2.getCrystallographicInfo();
						PDBHeader header2 = s2.getPDBHeader();
						SpaceGroup sp2 = info2.getSpaceGroup();
						String t2 = header2.getTechnique();
						float r2 = header2.getResolution();
						String c2 = header2.getClassification();
		//				System.out.println("Classification: " + header2.getClassification());
		//				System.out.println("Description: " + header2.getDescription());
		//				System.out.println("Title: " + header2.getTitle());
		//				AligQuaternaryStructure aligner = new AligQuaternaryStructure();
				 //       double[] rmsds = aligner.align(pdbIdI, pdbIdJ);
	//			        double[] rmsds = aligner.align(s1, s2);
	//			        double deltaRmsd = rmsds[0]-rmsds[1];
				        out.println(pdbId1 + "," + pdbId2 +"," + ba1 + "," + ba2 + "," + getFormula(structureIds.get(i)) + "," + 
	                            getPointGroup(structureIds.get(i)) + "," + 
				        		t1 + "," + t2 + "," + sp1 + "," + sp2 + "," + r1 + "," + r2 + "," + c1 + "," + c2);

				        out.flush();
						count++;
					}
				}
			}

		}
		out.close();
		error.close();
	}
	
	private String getBioassemblyId(String pdbId) {
		return csv.get(pdbId)[BIOASSEMBLY];
	}
	
	private String getFormula(String pdbId) {
		return csv.get(pdbId)[FORMULA];
	}
	
	private String getStoichiometry95(String pdbId) {
		return csv.get(pdbId)[STOICHIOMETRY95];
	}
	
	private String getPointGroup(String pdbId) {
		return csv.get(pdbId)[POINT_GROUP];
	}
	
	private void readFile(String fileName) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] tokens = line.split(",");
			String pdbId = tokens[PDB_ID];
			String bioassemblyId = tokens[BIOASSEMBLY];
			csv.put(pdbId+bioassemblyId, tokens);
		}
		reader.close();
	}
	
	private void createStructureGroups() {
		for (Entry<String, String[]> entry: csv.entrySet()) {
			String structureKey = entry.getKey();
			System.out.println("PDB " + structureKey + ": " + Arrays.toString(entry.getValue()));
			String[] tokens = entry.getValue();
			System.out.println("tokens: " + tokens.length);
			String formula = tokens[FORMULA];
			// exclude monomers for now
			if (formula.equals("A")) {
				continue;
			}
			System.out.println("formula: "+ formula);
			String stoichiometry95 = tokens[STOICHIOMETRY95];
			System.out.println("s95: " + stoichiometry95);
			if (! formula.equals(stoichiometry95)) {
				System.err.println(structureKey + ": inconsistent stoichiometry: " + formula + " - " + stoichiometry95);
				continue;
			}
			String signature95 = tokens[SIGNATURE95];
			String pointGroup = tokens[POINT_GROUP];
			String key = signature95 + "_" + pointGroup;
			System.out.println("Key: " + key);
		    if (structureGroups.containsKey(key)) {
		    	List<String> pdbIds = structureGroups.get(key);
		    	pdbIds.add(structureKey);
		    } else {
		    	List<String> list = new ArrayList<String>();
		    	list.add(structureKey);
		    	structureGroups.put(key, list);
		    }
		}
	}
	
	private Structure getStructure(String pdbId, PrintWriter error) {
		AtomCache cache = new AtomCache();
		cache.setAutoFetch(true);
		cache.setPath(PDB_PATH);	
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		p.setParseCAOnly(true);
		//p.setMaxAtoms(50000000);
		
		p.setAtomCaThreshold(Integer.MAX_VALUE);
	//	p.setAcceptedAtomNames(new String[]{" CA ", " CB "});
		cache.setFileParsingParams(p);

		Structure structure = null;
		try {
		    structure = cache.getStructure(pdbId);
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
