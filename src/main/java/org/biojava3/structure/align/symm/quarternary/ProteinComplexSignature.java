package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Chain;

public class ProteinComplexSignature {
	private BlastClustReader blastClust = null;
	private String pdbId = "";
	private List<Chain> chains = new ArrayList<Chain>();
	private List<List<Integer>> clusters = new ArrayList<List<Integer>>();
	private List<ChainSignature> chainSignatures = new ArrayList<ChainSignature>();
	

	public ProteinComplexSignature(String pdbId, GlobalSequenceGrouper grouper, BlastClustReader blastClust) {
		this.pdbId = pdbId;
		this.chains = grouper.getChains();
		this.clusters = grouper.getSequenceCluster100();
		this.blastClust = blastClust;
		
		createChainSignatures();
	}
	
	public String getComplexSignature() {
		StringBuilder builder = new StringBuilder();
		for (ChainSignature s: chainSignatures) {
			builder.append(s.toString());
		}
		return builder.toString();
	}
	
	public String getCompositionId(String chainId) {
		for (ChainSignature s: chainSignatures) {
			if (s.getChainIds().contains(chainId)) {
				return s.getCompositionId();
			}
		}
		return "";
	}
	
	private List<ChainSignature> createChainSignatures() {	
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	
		for (List<Integer> cluster: clusters) {
			String representativeChain = getRepresentativeChain(cluster.get(0));
			List<String> chainIds = new ArrayList<String>();
			for (Integer c: cluster) {
				chainIds.add(chains.get(c).getChainID());
			}
			
			ChainSignature chainSignature = new ChainSignature(representativeChain, cluster.size(), chainIds);
			String c = "?";
			if (clusters.indexOf(cluster) < alpha.length()) {
				c = alpha.substring(clusters.indexOf(cluster), clusters.indexOf(cluster)+1);
			}
			chainSignature.setCompositionId(c);
			System.out.println("Chains: " + chainIds + " compositionId: " + c);
            chainSignatures.add(chainSignature);
		}
		Collections.sort(chainSignatures);
		
		return chainSignatures;
	}

	
	
    private String getRepresentativeChain(int index) {
    	String chainId = chains.get(index).getChainID();
    	return blastClust.getRepresentativeChain(pdbId, chainId);
	}

}
