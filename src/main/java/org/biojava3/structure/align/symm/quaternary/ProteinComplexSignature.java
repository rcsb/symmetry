package org.biojava3.structure.align.symm.quaternary;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.biojava.bio.structure.Chain;

public class ProteinComplexSignature {
	private BlastClustReader blastClust = null;
	private String pdbId = "";
	private SequenceClusterer grouper = null;
	private List<Chain> chains = new ArrayList<Chain>();
	private List<List<Integer>> clusters = new ArrayList<List<Integer>>();
	private List<ChainSignature> chainSignatures = new ArrayList<ChainSignature>();
	

	public ProteinComplexSignature(String pdbId, SequenceClusterer grouper, BlastClustReader blastClust) {
		this.pdbId = pdbId;
		this.grouper = grouper;
		this.chains = grouper.getChains();
		this.clusters = grouper.getSequenceCluster100();
		this.blastClust = blastClust;
		
		getChainSignatures();
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
	
	public String getComplexStoichiometry() {
		StringBuilder s = new StringBuilder();
		for (ChainSignature c: chainSignatures) {
			s.append(c.getCompositionId());
			if (c.getChainIds().size() >1) {
				s.append(c.getChainIds().size());
			}
		}
		return s.toString();
	}
	
	public int getSubunitTypeCount() {
		return chainSignatures.size();
	}
	
	private List<ChainSignature> getChainSignatures() {	
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	
		Map<String,Integer> mapCounts = new TreeMap<String,Integer>();
		Map<String,List<String>> mapChainIds = new TreeMap<String, List<String>>();
		for (String chainId: grouper.getOrderedChainIDList()) {
			String rep = blastClust.getRepresentativeChain(pdbId, chainId);
			Integer value = mapCounts.get(rep);
			if (value == null) {
				mapCounts.put(rep, 1);
				List<String> list = new ArrayList<String>();
				list.add(chainId);
				mapChainIds.put(rep, list);
			} else {
				value+=1;
				mapCounts.put(rep, value);
				List<String> list = mapChainIds.get(rep);
				list.add(chainId);
			}
		}


		for (Entry<String, Integer> entry: mapCounts.entrySet()) {
			List<String> chainIds = mapChainIds.get(entry.getKey());
			ChainSignature chainSignature = new ChainSignature(entry.getKey(), entry.getValue(), chainIds);
			chainSignatures.add(chainSignature);
		}

		Collections.sort(chainSignatures);
		for (int i = 0; i < chainSignatures.size(); i++) {
			ChainSignature c = chainSignatures.get(i);
			if (i < alpha.length()) {
				c.setCompositionId(alpha.substring(i,i+1));
			} else {
				c.setCompositionId("?");
			}
		}
		
		return chainSignatures;
	}
	
	private List<ChainSignature> createChainSignaturesOld() {	
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
