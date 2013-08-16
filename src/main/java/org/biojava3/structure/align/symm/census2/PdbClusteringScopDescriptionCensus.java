/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2013-02-18
 *
 */
package org.biojava3.structure.align.symm.census2;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.MissingResourceException;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.align.client.FarmJobParameters;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.representatives.ScopSupport;

/**
 * A census that takes a list of sun Ids and a sequence clustering. Uses the sequence-clusteirng provided by RCSB.
 * 
 * @author dmyerstu
 */
public class PdbClusteringScopDescriptionCensus extends ScopDescriptionCensus {

	private static final Logger logger = LogManager.getLogger(PdbClusteringScopDescriptionCensus.class.getPackage().getName());

	private static final String SERVER_LOCATION = FarmJobParameters.DEFAULT_SERVER_URL;

	private int identityCutoff;

	public static void buildDefault(File censusFile, int identityCutoff, int[] sunIds) {
		try {
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			PdbClusteringScopDescriptionCensus census = new PdbClusteringScopDescriptionCensus(maxThreads, identityCutoff, sunIds);
			census.setOutputWriter(censusFile);
			census.setCache(new AtomCache());
			census.run();
			System.out.println(census);
		} catch (RuntimeException e) {
			logger.fatal(e.getMessage(), e);
		}
	}

	public static Set<String> getClusterRepresentatives(int cutoff) {
		Set<String> names = JFatCatClient.getRepresentatives(SERVER_LOCATION, cutoff);
		if (names == null) throw new MissingResourceException("Could not retrieve representatives from "
				+ SERVER_LOCATION, "JFatCatClient", "representative chains");
		return names;
	}

	public static void main(String[] args) {
		if (args.length < 2) {
			System.err.println("Usage: " + PdbClusteringScopDescriptionCensus.class.getSimpleName() + " census-file identity-cutoff sun-id-1 [sun-id-2 ...]");
			return;
		}
		ScopFactory.setScopDatabase(ScopFactory.getSCOP(ScopFactory.VERSION_1_75A));
		final File censusFile = new File(args[0]);
		final int identityCutoff = Integer.parseInt(args[1]);
		int[] sunIds = new int[args.length - 2];
		StringBuilder sb = new StringBuilder();
		for (int i = 2; i < args.length; i++) {
			Integer sunId = ScopSupport.getInstance().getSunId(args[i]);
			if (sunId == null) throw new IllegalArgumentException("Couldn't find " + args[i]);
			sb.append(sunId + ",");
			sunIds[i - 2] = sunId;
		}
		logger.info("Running on " + sb.toString().substring(0, sb.toString().length()-1));
		buildDefault(censusFile, identityCutoff, sunIds);
	}

	/**
	 * We want to match domains to particular chains Unfortunately, chains can contain more than one domain, and domains
	 * can be defined over more than one chain. Since for any chains A and B for some PDB id we won't have both a domain
	 * defined over A AND a domain defined over A and B (at least I hope not), it's safe to include any domain that is
	 * defined over a chain that is a cluster representative, even if the domain also includes other chains.
	 */
	public static boolean isDomainOverChain(ScopDomain domain, Set<String> clusterRepresentatives) {
		final List<String> ranges = domain.getRanges();
		for (String range : ranges) {
			final int index = range.indexOf(':');
			final String chain;
			if (index == -1) {
				// this happens when the string is just - (meaning everything/all?)
				// ordinarily this fallback wouldn't work in every case
				// but it's okay when we have -
				final StructureName scopName = new StructureName(domain.getScopId());
				chain = scopName.getPdbId() + "." + scopName.getChainId();
			} else {
				chain = domain.getPdbId() + "." + range.substring(0, index);
			}
			// System.out.println(chain);
			if (clusterRepresentatives.contains(chain.toUpperCase())) { // the toUpperCase is critical
				return true;
			}
		}
		return false;
	}

	public PdbClusteringScopDescriptionCensus(int maxThreads, int identityCutoff, int[] sunIds) {
		super(maxThreads, sunIds);
		this.identityCutoff = identityCutoff;
	}

	@Override
	protected List<ScopDomain> getDomains() {

		final Set<String> clusterRepresentatives = getClusterRepresentatives(identityCutoff);
		List<ScopDomain> domains = new ArrayList<ScopDomain>();

		for (int sunId : sunIds) {

			List<ScopDomain> putative = new ArrayList<ScopDomain>();
			ScopDescriptionCensus.getDomainsUnder(sunId, putative);

			for (ScopDomain domain : putative) {
				if (isDomainOverChain(domain, clusterRepresentatives)) domains.add(domain);
			}

		}

		return domains;
	}

}
