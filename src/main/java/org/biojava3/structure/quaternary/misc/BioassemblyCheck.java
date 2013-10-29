package org.biojava3.structure.quaternary.misc;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.biojava3.structure.quaternary.utils.BlastClustReader;

public class BioassemblyCheck {
	private static String RESULT_DIR = "C:/Users/Peter/Documents/QuatStructureComparison/";
	
	private List<String> pdbids = Collections.emptyList();
	private List<String> signatures =  Collections.emptyList();
	private List<String> stoichiometries =  Collections.emptyList();
	private List<String> pointGroups =  Collections.emptyList();
	private List<String> localSymmetry = Collections.emptyList();
	private List<String> pseudoSymmetry = Collections.emptyList();
	
	private static int sequenceIdentity = 95;
	private static String PDB_PATH = "C:/PDB/";
	private static String fileName = "C:/Users/Peter/Desktop/Composition Symmetry/rep_sym 90 20120831.csv";
	private static String outFile = "C:/Users/Peter/Desktop/Composition Symmetry/bioassemblyCheck";

	private BlastClustReader reader95 = null;
	
	private Map<String, List<Integer>> entityMap = new HashMap<String, List<Integer>>();
	private Map<String, Integer> entityFormulaMap = new TreeMap<String, Integer>();
	private Map<String, Integer> entityFormulaPointgroupMap = new TreeMap<String, Integer>();
	private List<List<String>> table = new ArrayList<List<String>>();
	
	public BioassemblyCheck() {
		reader95 = new BlastClustReader(95);
		SimpleCsvReader reader = new SimpleCsvReader();
		try {
			reader.readFile(RESULT_DIR + "20131007_205827_symm.csv");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		pdbids = reader.getColumn("pdbId");
		signatures = reader.getColumn("signature95");
		stoichiometries = reader.getColumn("stoichiometry");
		pointGroups = reader.getColumn("pointgroup");
		localSymmetry = reader.getColumn("local");
		pseudoSymmetry = reader.getColumn("pseudosymmetric");
		createEntityMap();
	}
	
	private void createEntityMap() {

		for (int i = 0; i < pdbids.size(); i++) {
			String entityList = getEntities(signatures.get(i));
			List<Integer> indices = entityMap.get(entityList);
			if (indices == null) {
				indices = new ArrayList<Integer>();
				entityMap.put(entityList, indices);
			}
			indices.add(i);
		}
	}
	
	public void BioassemblyCheck(List<PdbBlastHit> blastHits, String stoichiometry, String pointGroup) {
		System.out.println("BlastHits: " + blastHits);
	
		String signature = getSignature(blastHits);
		System.out.println("Blast hit signature: " + signature);

		Map<String, List<String>> countMap = new HashMap<String, List<String>>();
		if (entityMap.containsKey(signature)) { 
			List<Integer> indices = entityMap.get(signature);
			for (int index: indices) {
				if (localSymmetry.get(index).equals("false") && pseudoSymmetry.get(index).equals("false")) {
					String key = stoichiometries.get(index) + " " + pointGroups.get(index);
					List<String> pdbList = countMap.get(key);
					if (pdbList == null) {
						pdbList = new ArrayList<String>();
						countMap.put(key, pdbList);
					} 
					pdbList.add(pdbids.get(index).substring(3));
				}
			}
		}
		
		System.out.println();
		System.out.println("Stoichiometry and Point Group statistics: ");
		int total = 0;
		for (Entry<String, List<String>> entry: countMap.entrySet()) {
			total+= entry.getValue().size();
			System.out.println(entry.getKey() + "(" + entry.getValue().size() + "): " + entry.getValue());
		}
		System.out.println();
		String sampleKey = stoichiometry + " " + pointGroup;
		List<String> pdbList = countMap.get(sampleKey);
		int count = pdbList == null ? 0 : pdbList.size();
		float score = (float)count/total;
		System.out.println("Score for sample " + sampleKey + ": " + score + " (" + count + "/" + total + ")");
		System.out.println();
		System.out.println();
	}

	private String getSignature(List<PdbBlastHit> blastHits) {
		List<String> signatures = new ArrayList<String>();
		
		for (PdbBlastHit hit: blastHits) {
			String rep = reader95.getRepresentativeChain(hit.getPdbId(), hit.getChainIds().get(0));
			System.out.println("Representative for: " + hit.getPdbId()+hit.getChainIds().get(0) + " is " + rep);
			signatures.add(rep);
		}
		Collections.sort(signatures);
		String signature = "";
		for (int i = 0, n = signatures.size(); i < n; i++) {
			signature += signatures.get(i);
			if (i < n-1) {
				signature += "_";
			}
		}
		return signature;
	}
	
	
	
//	private void run() throws IOException {
//		readTable(fileName);
//		calcMap(entityFormulaMap, 0);
//		calcMap(entityFormulaPointgroupMap, 1);
//		calcMap(entityMap, 2);
//		printMap();
//		calcStats();
//		printList();
//		saveResults(outFile);
//	}
	
	private void calcMap(Map<String, Integer> map, int index) {
		String key = "";
		for (List<String> row: table) {
			key = row.get(index);
			Integer count = map.get(key);
			if (count == null) {
				count = new Integer(0);
			}
			count++;
			
			map.put(key, count);
		}
	}

	private void calcStats() {
//		for (List<String> row: table) {
//			String key = row.get(0);
//			Integer subTotalFormula = entityFormulaMap.get(key);
//			key = row.get(1);
//			Integer subTotalFormulaPointgroup = entityFormulaPointgroupMap.get(key);
//			key = row.get(2);
//			Integer total = entityMap.get(key);
//			float ratio = (float)subTotalFormula/(float)total;
//			row.add(Float.toString(ratio));
//			ratio = (float)subTotalFormulaPointgroup/(float)total;
//			row.add(Float.toString(ratio));
//		}
	}

	private void printMap() {
		for (Entry<String, Integer> iter: entityFormulaMap.entrySet()) {
				System.out.println(iter.getKey() + ": " + iter.getValue());
		}
	}
	
	private void printList() {
		for (List<String> row: table) {
				System.out.println(row);
		}
	}
	
	private static String getEntities(String signature) {
		StringBuilder b = new StringBuilder();
		int from = 0;
		do {
			from = signature.indexOf("(", from);
			if (from < 0) break;
			if (from > 0) b.append("_");
			int end = signature.indexOf(")", from+1);
			b.append(signature.substring(from+1, end)) ;
			from = end + 1;
		} while (true);
		return b.toString();
	}
	
	private void saveResults(String outFile) throws IOException {
		PrintWriter out = new PrintWriter(outFile+ sequenceIdentity + ".csv");
		out.println("Sequence cluster " + sequenceIdentity + "%,PDB ID,Bioassemby,Formula,Point group,Consistency,Consistency_sym");
		for (List<String> row: table) {
			row = row.subList(2, row.size());
			if (! (row.get(5).equals("1.0") && row.get(6).equals("1.0"))) {
				String s = row.toString();
				out.println(row.toString().substring(1,s.length()-1));
			}
		}
		out.flush();
		out.close();
	}
}
