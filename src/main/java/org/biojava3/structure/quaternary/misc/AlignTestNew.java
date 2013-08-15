package org.biojava3.structure.quaternary.misc;

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
import org.biojava3.structure.quaternary.core.ChainClusterer;
import org.biojava3.structure.quaternary.core.ClusterProteinChains;
import org.biojava3.structure.quaternary.core.ProteinSequenceClusterer;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.SequenceAlignmentCluster;
import org.biojava3.structure.quaternary.core.Subunits;

public class AlignTestNew {

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
	
		
		AlignTestNew test = new AlignTestNew();
		test.cluster("4HHB","3R5I");

	}
	
	public void cluster(String pdbId1, String pdbId2) {
		Structure s1 = getStructure(pdbId1);
		Structure s2 = getStructure(pdbId2);
		QuatSymmetryParameters parameters = new QuatSymmetryParameters();
		
		ClusterProteinChains clusterer = new ClusterProteinChains(s1, s2, parameters);
		ChainClusterer chainClusterer = new ChainClusterer(clusterer.getSequenceAlignmentClusters(0.95));
		System.out.println(chainClusterer);
		
		Subunits subunits = new Subunits(chainClusterer.getCalphaCoordinates(), 
				chainClusterer.getSequenceClusterIds(),
				chainClusterer.getPseudoStoichiometry(),
				chainClusterer.getMinSequenceIdentity(),
				chainClusterer.getMaxSequenceIdentity(),
				chainClusterer.getFolds(),
				chainClusterer.getChainIds(),
				chainClusterer.getModelNumbers());
	}
		
	private Structure getStructure(String pdbId) {
		AtomCache cache = new AtomCache();
		cache.setAutoFetch(true);	
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);	
		p.setAtomCaThreshold(Integer.MAX_VALUE);
		cache.setFileParsingParams(p);
			
		Structure structure = null;
		try {
			structure = cache.getBiologicalAssembly(pdbId, 1, true);
		} catch (StructureException e) {
			e.printStackTrace();
			System.err.println(pdbId + "------------------------------------");
			System.err.println(e.getMessage());
			System.err.flush();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println(pdbId + "------------------------------------");
			System.err.println(e.getMessage());
			System.err.flush();
		}
		return structure;
	}
}
