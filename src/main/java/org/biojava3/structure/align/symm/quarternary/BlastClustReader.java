package org.biojava3.structure.align.symm.quarternary;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.align.xml.RepresentativeXMLConverter;
import org.rcsb.fatcat.server.PdbChainKey;

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
			URLConnection connection = u.openConnection();
			connection.setConnectTimeout(60000);
			InputStream stream = connection.getInputStream();

			if (stream != null) {
				BufferedReader reader = new BufferedReader(new InputStreamReader(stream));

				String line = null;
				try {
					while ((line = reader.readLine()) != null) {
						line = line.replace('_', '.');
						List<String> cluster = Arrays.asList(line.split(" "));	
						clusters.add(cluster);
					}
				} catch (IOException e) {
					//e.printStackTrace();
				} finally {
					try {
						stream.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		return;
	}
	
}

