package org.biojava3.structure.align.symm.quaternary.analysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

public class BioassemblyCheck {
	private static int sequenceIdentity = 95;
	private static String PDB_PATH = "C:/PDB/";
	private static String fileName = "C:/Users/Peter/Desktop/Composition Symmetry/rep_sym 90 20120831.csv";
	private static String outFile = "C:/Users/Peter/Desktop/Composition Symmetry/bioassemblyCheck";
	private int pdbIdIndex = -1;
	private int bioassemblyIndex = -1;
	private int formulaIndex = -1;
	private int signatureIndex = -1;
	private int pointgroupIndex = -1;
	
	private Map<String, Integer> entityMap = new TreeMap<String, Integer>();
	private Map<String, Integer> entityFormulaMap = new TreeMap<String, Integer>();
	private Map<String, Integer> entityFormulaPointgroupMap = new TreeMap<String, Integer>();
	private List<List<String>> table = new ArrayList<List<String>>();
	
	public static void main(String[] args) throws IOException {
		new BioassemblyCheck().run();
	}
	
	private void run() throws IOException {
		readTable(fileName);
		calcMap(entityFormulaMap, 0);
		calcMap(entityFormulaPointgroupMap, 1);
		calcMap(entityMap, 2);
		printMap();
		calcStats();
		printList();
		saveResults(outFile);
	}

	private void readTable(String fileName) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));

		String line = reader.readLine();

		// find index of items in header
		List<String> items = Arrays.asList(line.split(","));	
		System.out.println(items);
		for (int i = 0; i < items.size(); i++) {
			if (items.get(i).equals("pdbId")) {
				pdbIdIndex = i;
			}
			if (items.get(i).equals("bioassembly")) {
				bioassemblyIndex = i;
			}
			if (items.get(i).equals("formula")) {
				formulaIndex = i;
			}
			if (items.get(i).equals("signature"+sequenceIdentity)) {
				signatureIndex = i;
			}
			if (items.get(i).equals("pointgroup")) {
				pointgroupIndex = i;
			}
		}

		while ((line = reader.readLine()) != null) {
			//			System.out.println(items);
			items = Arrays.asList(line.split(","));
			String pdbId = items.get(pdbIdIndex);
			String bioassembly = items.get(bioassemblyIndex);
			String formula = items.get(formulaIndex);
			String entities = getEntities(items.get(signatureIndex));
			String pointgroup = items.get(pointgroupIndex);
			String key1 = entities + "_" + formula;
			String key2 = entities + "_" + formula + "_" + pointgroup;
			List<String> list = new ArrayList<String>();
			list.add(key1);
			list.add(key2);
			list.add(entities);
			list.add(pdbId);
			list.add(bioassembly);
			list.add(formula);
			list.add(pointgroup);
	
			table.add(list);
		}
	}
	
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
		for (List<String> row: table) {
			String key = row.get(0);
			Integer subTotalFormula = entityFormulaMap.get(key);
			key = row.get(1);
			Integer subTotalFormulaPointgroup = entityFormulaPointgroupMap.get(key);
			key = row.get(2);
			Integer total = entityMap.get(key);
			float ratio = (float)subTotalFormula/(float)total;
			row.add(Float.toString(ratio));
			ratio = (float)subTotalFormulaPointgroup/(float)total;
			row.add(Float.toString(ratio));
		}
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
