package org.biojava3.structure.dbscan;

import java.io.InputStream;

import java.net.URL;

import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;

import java.util.TreeSet;

import org.biojava.bio.structure.align.client.JFatCatClient;

import org.biojava.bio.structure.align.util.HTTPConnectionTools;

import org.biojava.bio.structure.align.xml.RepresentativeXMLConverter;

import org.rcsb.fatcat.server.PdbChainKey;

public class GetRepresentatives {

	private static String coreUrl = "http://www.rcsb.org/pdb/rest/representatives?cluster=";
	// available sequence clusters
	private static List<Integer> seqIdentities = Arrays.asList(30, 40, 50, 70, 90, 95, 100);

	/**
	 * Returns a representative set of PDB protein chains at 40% sequence 
	 * identity cutoff.
	 * @param sequenceIdentity
	 * @return PdbChainKey set of representatives
	 * @deprecated
	 */
	public static SortedSet<PdbChainKey> getRepresentatives() {
        return getRepresentatives(40);
	}

	/**
	 * Returns a representative set of PDB protein chains at the specified sequence 
	 * identity cutoff. See http://www.pdb.org/pdb/statistics/clusterStatistics.do
	 * for more information.
	 * @param sequenceIdentity sequence identity threshold
	 * @return PdbChainKey set of representatives
	 */
	public static SortedSet<PdbChainKey> getRepresentatives(int sequenceIdentity) {
		SortedSet<PdbChainKey> representatives = new TreeSet<PdbChainKey>();

		if (!seqIdentities.contains(sequenceIdentity)) {
			System.err.println("Error: representative chains are not available for %sequence identity: "
							+ sequenceIdentity);
			return representatives;
		}


		try {

			URL u = new URL(coreUrl + sequenceIdentity);

			InputStream stream = HTTPConnectionTools.getInputStream(u, 60000);

			String xml = null;

			if (stream != null) {
				xml = JFatCatClient.convertStreamToString(stream);

				SortedSet<String> reps = RepresentativeXMLConverter.fromXML(xml);

				for (String s : reps) {
					PdbChainKey k = PdbChainKey.fromName(s);
					representatives.add(k);
				}

			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		return representatives;
	}

}
