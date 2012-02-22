package org.biojava3.structure.align.symm.quarternary;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;

public class AligQuaternaryStructure {
	private static String PDB_PATH = "C:/PDB/";
	private static final int MIN_SEQUENCE_LENGTH = 24;
	private List<String> pdbIds = new ArrayList<String>();
	private Map<Integer,Integer> clusterMap = new LinkedHashMap<Integer,Integer>();
	private List<Atom[]> ca1 = null;
	private List<Atom[]> ca2 = null;
	private List<Atom> caQuaternary1 = new ArrayList<Atom>();
	private List<Atom> caQuaternary2 = new ArrayList<Atom>();

	public AligQuaternaryStructure(List<String> pdbIds) {
		this.pdbIds = pdbIds;
	}
	
	public void align() {
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
		
		PrintWriter out = null;
		PrintWriter error = null;
		try {
			out = new PrintWriter(new FileWriter(PDB_PATH + "quaternaryAlignment"));
			error = new PrintWriter(new FileWriter(PDB_PATH + "error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		for (int i = 0; i < pdbIds.size(); i++) {
			String pdbId1 = pdbIds.get(i);
			Structure s1 = getStructure(pdbId1);
			for (int j = 0; j < pdbIds.size(); j++) {
				String pdbId2 = pdbIds.get(j);
				Structure s2 = getStructure(pdbId2);
				System.out.println("----------- Aligning: " + pdbId1 + " - " + pdbId2);
				alignUniqueSequences(s1, s2);
			}
		}

//		Structure s1 = getStructure(pdbId1);
//		System.out.println(s1.getChains().size());
//		Structure s2 = getStructure(pdbId2);
//		System.out.println(s2.getChains().size());
//
//		alignUniqueSequences(s1, s2);
	}
	
	private void alignUniqueSequences(Structure s1, Structure s2) {
		GlobalSequenceGrouper grouper1 = new GlobalSequenceGrouper(s1, MIN_SEQUENCE_LENGTH);
		List<List<Integer>> clusters1 = grouper1.getSequenceCluster100();
        ca1 = grouper1.getCalphaTraces();
        System.out.println("Cluster1: " + clusters1.size());

        GlobalSequenceGrouper grouper2 = new GlobalSequenceGrouper(s2, MIN_SEQUENCE_LENGTH);
        List<List<Integer>> clusters2 = grouper2.getSequenceCluster100();
        ca2 = grouper2.getCalphaTraces();
        System.out.println("Cluster2: " + clusters1.size());
        clusterMap = new LinkedHashMap<Integer,Integer>();
        
        int[][][] alig = null;
        
        for (List<Integer> c1: clusters1) {
        	// get the first C alpha array from a cluster, as the representative of this cluster
        	int representative1 = c1.get(0); 

        	Atom[] ca1Seq = ca1.get(representative1);
        	for (List<Integer> c2: clusters2) {
        		int representative2 = c2.get(0); 
            	System.out.print("Comparing: " + representative1);
            	System.out.println(" - " + representative2);
            	Atom[] ca2Seq = ca2.get(representative2);
        		double identity = alignPair(ca1Seq, ca2Seq, alig);
        		if (identity > 0.99) {
        			clusterMap.put(representative1,representative2);
        			createCalphaList(ca1Seq, ca2Seq, alig);
        		}
        	}
        }
 
        for (Entry<Integer, Integer> e: clusterMap.entrySet()) {
        	System.out.println(e);
        }
	}
	
	private void createCalphaList(Atom[] ca1Seq, Atom[] ca2Seq, int[][][] alig) {
		int[] numbering1 = alig[0][0];
		int[] numbering2 = alig[0][1];
		for (int i = 0; i < numbering1.length; i++) {
			
		}
		
	}

	private double calcQuaternaryStructureRmsd() {
		double rmsd = 0.0;
		
		return rmsd;
	}
	
	private double alignPair(Atom[] ca1Seq, Atom[] ca2Seq, int[][][] alig) {
		SmithWaterman3Daligner aligner = new SmithWaterman3Daligner();
		AFPChain afp = null;
		try {
			afp = aligner.align(ca1Seq, ca2Seq);
		} catch (StructureException e) {
			e.printStackTrace();
			return 0.0f;
		} 
		alig = afp.getOptAln();

		double identity = afp.getIdentity();
		if (identity > 0.99) {
			System.out.println("Aligned residues: ");
			System.out.println(Arrays.toString(alig[0][0]));
			System.out.println(Arrays.toString(alig[0][1]));
			System.out.println("Aligned residues: " + alig[0][0].length);
			System.out.println("Seq. identity: " + afp.getIdentity());
			System.out.println("Rmsd: " + afp.getChainRmsd());
		}
		return identity;
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
