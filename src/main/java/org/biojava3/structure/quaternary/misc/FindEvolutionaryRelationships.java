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
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class FindEvolutionaryRelationships {
	private static String PDB_PATH = "C:/PDB/";
	private static String fileName = "c:/Users/Peter/Desktop/rep_sym 20120317.csv";
	private Map<String, String[]> csv = new HashMap<String,String[]>();
	private List<String> uniqueSignatures = new ArrayList<String>();
	private List<String> stoichiometries = new ArrayList<String>();
	private Map<String,List<String>> chainToSignatures = new TreeMap<String,List<String>>();
	private List<Map<String, Integer>> compositions = new ArrayList<Map<String, Integer>>();
	private List<List<String>> relationships = new ArrayList<List<String>>();
	private int sequenceIdentity = 100;

	public static void main(String[] args) {
		FindEvolutionaryRelationships finder = new FindEvolutionaryRelationships();
		finder.run(100);
	    finder = new FindEvolutionaryRelationships();
		finder.run(90);
		finder = new FindEvolutionaryRelationships();
		finder.run(70);
		finder = new FindEvolutionaryRelationships();
		finder.run(40);
	}
	
	public void run(int sequenceIdentity) {
		try {
			readFile(fileName);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.sequenceIdentity = sequenceIdentity;
		uniqueSignatures(sequenceIdentity);
		for (String s: uniqueSignatures) {
//			System.out.println(s);
		}
        mapSignaturesToChains();
		analyze();
		aggregateResults("quatNet" + sequenceIdentity);
		System.out.println("Unique signatures: (" + sequenceIdentity + "): " + uniqueSignatures.size());
	}
	
	private void test() {
//		String signature1 = "(1A03_A)(2JTT_B)";
//		String signature2 = "(2JTT_B)2(1A03_A)2";
//		String signature1 = "(1P5U_A)(1P5U_B)(1P5U_B)";
//		String signature2 = "(1P5U_A)(1P5U_B)";
		String signature1 = "(1VYW_A)(1VYW_B)";
		String signature2 = "(1F5Q_B)(1VYW_A)";
//		System.out.println("isSubstucture: " + isSubstructure(signature1, signature2));
		System.out.println("Unique signatures: " + uniqueSignatures.size());
	}
	private void uniqueSignatures(int sequenceIdentity) {
		int tmp = 0;
			for (Entry<String,String[]> entry: csv.entrySet()) {
			String signature = getSignature(entry, sequenceIdentity);
			// skip incomplete signatures
			if (signature.startsWith("()")) {
				continue;
			}
			// skip header
			if (signature.startsWith("signature")) {
				continue;
			}
			if (!uniqueSignatures.contains(signature)) {
				// composition map will be empty for entries that have duplicate chain ids,
				// this happens when there is a mismatch between the uniprot seq. and the atom seq.,
				// i.e., when modified residues are in the sequence
				Map<String,Integer> composition = getCompositionMap(signature);
				if (!composition.isEmpty()) {
					compositions.add(composition);
					uniqueSignatures.add(signature);
					stoichiometries.add(getStoichiometry(entry, sequenceIdentity));
				}
			}
			tmp++;
//			if (tmp > 500) {
//				return;
//			}
		}
	}
	
	private void mapSignaturesToChains() {
		for (String signature: uniqueSignatures) {	
			Map<String,Integer> composition = getCompositionMap(signature);
			for (String chainId: composition.keySet()) {
				List<String> s = chainToSignatures.get(chainId);
				if (s == null) {
					s = new ArrayList<String>();
					chainToSignatures.put(chainId, s);
	
				}
				s.add(signature);
			}
		}
	}
	
	private Map<String,Integer> getCompositionMap(String signature) {
//		System.out.println("sign: " + signature);
		String components = signature;
		components = signature.replaceAll("\\(", " ");
		components = components.replaceAll("\\)", " ");
		components = components.replaceAll("  ", " ");;
//		System.out.println("split: " + components);
		String[] parts = components.split(" ");
//		System.out.println("decomp: "  + Arrays.toString(parts));
		Map<String, Integer> compositionMap = new HashMap<String, Integer>();
		String chainId = "";
		int multiplier = 1;
		for (int i = 0; i < parts.length; i++) {
			if (parts[i].isEmpty()) {
				continue;
			}
//			System.out.println(parts[i]);
			// can be length 6 or 7
			if (parts[i].length() >= 6) {
				chainId = parts[i];
			}
			
			if (i < parts.length-1 && parts[i+1].length() < 6) {
				i++;
				if (parts[i].isEmpty()) {
					continue;
				}
				multiplier = Integer.parseInt(parts[i]);
			} else {
				multiplier = 1;
			}
			if (compositionMap.containsKey(chainId)) {
				System.out.println("Duplicate chain ids in signature: " + signature);
				return new HashMap<String, Integer>();
			} else {
				compositionMap.put(chainId, multiplier);
			}
//			System.out.println("parts: " + chainId + " x " + multiplier);
		}
		return compositionMap;
	}
	
	private void analyze() {
		for (int i = 0; i < compositions.size(); i++) {
			for (int j = i + 1; j < compositions.size(); j++)   
				if (isSubstructure(i, j) == 1) {
					//        			System.out.println("Substructure: " + s1 + " -> " + s2);
					List<String> data = new ArrayList<String>();
					data.add(stoichiometries.get(i));
					data.add(stoichiometries.get(j));
					relationships.add(data);
	//				System.out.println(data);
				}
		}
	}
	
	private int isSubstructure(int i, int j) {
		Map<String,Integer> map1 = compositions.get(i);
		Map<String,Integer> map2 = compositions.get(j);
		
		if (map1.size() > map2.size()) {
			return 0; // can't be a substructure or exact match
		}
		
		int substructure = 0;
		int exact = 0;
		for (Entry<String,Integer> entry1: map1.entrySet()) {
			String key1 = entry1.getKey();
			int mul1 = entry1.getValue();
			if (!map2.containsKey(key1)) {
				return 0; // can't be a substructure or exact match
			}
			for (Entry<String,Integer> entry2: map2.entrySet()) {
				String key2 = entry2.getKey();
				int mul2 = entry2.getValue();
				if (key1.equals(key2)) {
				  if ( mul1 < mul2) {
					  substructure++;
				  } else if (mul1 == mul2) {
					  exact++;
				  }
				}
			}
		}
//		System.out.println("map1: " + map1.size());
//		System.out.println("map2: " + map2.size());
//		System.out.println("substructures: " + substructure);
//		System.out.println("exact: " + exact);
		
		if (exact == map2.size() && map1.size() == map2.size()) {
			return 2; // exact match
		}
	    if (substructure > 0 || exact > 0) {
			return 1; // substructure match
		} else {
			return 0;
		}
	}
	
	private void aggregateResults(String filename) {
		int thresholdCount = 20;
		Map<String,Integer> count = new LinkedHashMap<String,Integer>();
		Map<String,String> composition1 = new LinkedHashMap<String,String>();
		Map<String,String> composition2 = new LinkedHashMap<String,String>();
		
		for (List<String> list: relationships) {
			String c1 = list.get(0);
			String c2 = list.get(1);
			String key = c1 + "_" + c2;
			Integer value = count.get(key);
			if (value == null) {
				value = new Integer(1);
				count.put(key, value);
				composition1.put(key, c1);
				composition2.put(key, c2);
//				System.out.println("Adding: " + key + "," + c1 +"," + c2 + "," + value);
			} else {
				value++;
				count.put(key, value);
//				System.out.println("Incr:  " + key + "," + c1 +"," + c2 + "," + value);
			}
			
		}
		 
		PrintWriter out = null;
		PrintWriter out1 = null;
		try {
			out = new PrintWriter(new FileWriter(PDB_PATH + filename + "_" + thresholdCount + ".csv"));
	//		out1 = new PrintWriter(new FileWriter(PDB_PATH +"_homo_homo_" + filename + ".csv"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		out.println("Stoichiometry1,Stoichiometry2,count");
		for (Entry<String, Integer> entry: count.entrySet()) {
			String key = entry.getKey();
	//		System.out.println("key: " + key);
			int value = entry.getValue();
			String c1 = composition1.get(key);
			String c2 = composition2.get(key); 
			if (value >= thresholdCount) {
			out.println(c1 + "," + c2 + "," + value);
			// homo-oligomers only
			if (! (c1.contains("B") || c2.contains("B"))) {
		//	    out1.println(c1 + "," + c2 + "," + value);
			}
			}
		}
		out.flush();
		out.close();
		//out1.flush();
		//out1.close();
	}
	
	private String getComposition(String pdbId) {
		return csv.get(pdbId)[1];
	}
	
	private String getSignature(Entry<String,String[]> entry, int sequenceIdentity) {
		if (sequenceIdentity == 100) {
		    return entry.getValue()[2];
		} else if (sequenceIdentity == 90) {
				return entry.getValue()[5];
		} else if (sequenceIdentity == 70) {
			return entry.getValue()[8];
		} else if (sequenceIdentity == 40) {
		    return entry.getValue()[11];
	    }
		return "";
	}
	
	private String getStoichiometry(Entry<String,String[]> entry, int sequenceIdentity) {
		  return entry.getValue()[3]; // try using 100% id as a reference stoichiometry
//		if (sequenceIdentity == 100) {
//		    return entry.getValue()[3];
//		} else if (sequenceIdentity == 90) {
//				return entry.getValue()[6];
//		} else if (sequenceIdentity == 70) {
//			return entry.getValue()[9];
//		} else if (sequenceIdentity == 40) {
//		    return entry.getValue()[12];
//	    }
//		return "";
	}
	
	private String getPointGroup(Entry<String,String[]> entry) {
		return entry.getValue()[14];
	}
	private String getLigands(String pdbId) {
		String[] tokens = csv.get(pdbId);
		if (tokens.length == 18) {
			return tokens[17];
		} else {
			return "";
		}
	}
	
	private void readFile(String fileName) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] tokens = line.split(",");
			String pdbId = tokens[0];
			csv.put(pdbId, tokens);
		}
		reader.close();
	}
}
