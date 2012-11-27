package org.biojava3.structure.quaternary.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


public class BlastClustReader {
	private int sequenceIdentity = 0;
	private List<List<String>> clusters = new ArrayList<List<String>>();
	private static final String coreUrl = "ftp://resources.rcsb.org/sequence/clusters/";
	private static List<Integer> seqIdentities = Arrays.asList(30, 40, 50, 70, 90, 95, 100);

	public BlastClustReader(int sequenceIdentity) {
		this.sequenceIdentity = sequenceIdentity;
	}
	
	public List<List<String>> getPdbChainIdClusters() {
		loadClusters(sequenceIdentity);
		return clusters;
	}
	
	public Map<String,String> getRepresentatives(String pdbId) {
		loadClusters(sequenceIdentity);
		String pdbIdUc = pdbId.toUpperCase();
		
		Map<String,String> representatives = new LinkedHashMap<String,String>();
		for (List<String> cluster: clusters) {
			// map fist match to representative
			for (String chainId: cluster) {
				if (chainId.startsWith(pdbIdUc)) {
					representatives.put(chainId, cluster.get(0));
                    break;
				}
			}
		}
		return representatives;
	}
	
	public String getRepresentativeChain(String pdbId, String chainId) {
		loadClusters(sequenceIdentity);
	
		// check if chain id is lower case. In that case double the chain id, i.e. "o" becomes "oo" (example PDB 1VU1_oo)
		// This appears to be the convention in the BlastClust files for lower case letters
		String cId = new String(chainId);
        String chainIdLc = cId.toLowerCase();
//        if (chainIdLc.equals(chainId) && Character.isAlphabetic(chainId.codePointAt(0))) { // not available in Java 5
        if (chainIdLc.equals(chainId) && !Character.isDigit(chainId.codePointAt(0))) {
        	cId = pdbId + "_" + chainIdLc + chainIdLc;
        } else {
        	cId = pdbId + "_" + chainId;
        }
		
		for (List<String> cluster: clusters) {
			for (String chnId: cluster) {
				if (chnId.equals(cId)) {
					return cluster.get(0);
				}
			}
		}
		return "";
	}
	
	public List<List<String>> getPdbChainIdClusters(String pdbId) {
		loadClusters(sequenceIdentity);
		String pdbIdUc = pdbId.toUpperCase();

		List<List<String>> matches = new ArrayList<List<String>>();
		for (List<String> cluster: clusters) {
			for (String chainId: cluster) {
				if (chainId.startsWith(pdbIdUc)) {
					matches.add(cluster);
					break;
				}
			}
		}
		return matches;
	}
	
	public List<List<String>> getChainIdsInEntry(String pdbId) {
		loadClusters(sequenceIdentity);
		String pdbIdUc = pdbId.toUpperCase();
		
		List<List<String>> matches = new ArrayList<List<String>>();
		List<String> match = null;
		
		for (List<String> cluster: clusters) {
			for (String chainId: cluster) {
				if (chainId.startsWith(pdbIdUc)) {
					if (match == null) {
						match = new ArrayList<String>();
					}
					match.add(chainId.substring(5));
				}
			}
			if (match != null) {
				Collections.sort(match);
				matches.add(match);
				match = null;
			}
		}
		return matches;
	}
	
	private void loadClusters(int sequenceIdentity) {
		// load clusters only once
		if (clusters.size() > 0) {
			return;
		}

		if (!seqIdentities.contains(sequenceIdentity)) {
			System.err.println("Error: representative chains are not available for %sequence identity: "
					+ sequenceIdentity);
			return;
		}

		try {
			URL u = new URL(coreUrl + "bc-" + sequenceIdentity + ".out");
			InputStream stream = u.openStream();
	//		URLConnection connection = u.openConnection();
	//		connection.setConnectTimeout(60000);
	//		InputStream stream = connection.getInputStream();

			if (stream != null) {
				BufferedReader reader = new BufferedReader(new InputStreamReader(stream));

				String line = null;
				try {
					while ((line = reader.readLine()) != null) {
	//					line = line.replace('_', '.');
						List<String> cluster = Arrays.asList(line.split(" "));	
						clusters.add(cluster);
					}
					reader.close();
					stream.close();
				} catch (IOException e) {
					//e.printStackTrace();
				} finally {
//					try {
//						System.out.println("closing reader");
//						reader.close();
//						stream.close();
//					} catch (IOException e) {
//						e.printStackTrace();
//					}
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return;
	}
	
}

