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

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.client.FarmJobParameters;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * A census that takes a list of sun Ids and a sequence clustering.
 * 
 * @author dmyerstu
 */
public class AstralScopDescriptionCensus extends ScopDescriptionCensus {

	private static final String SERVER_LOCATION = FarmJobParameters.DEFAULT_SERVER_URL;

	static final Logger logger = Logger.getLogger(CensusJob.class.getPackage().getName());
	private int identityCutoff;
	private int[] sunIds;

	static {
		BasicConfigurator.configure();
	}

	public static void buildDefault(String pdbDir, File censusFile, int identityCutoff, int[] sunIds) {
		try {
			Census.setBerkeleyScop(pdbDir);
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			AstralScopDescriptionCensus census = new AstralScopDescriptionCensus(maxThreads, identityCutoff, sunIds);
			census.setOutputWriter(censusFile);
			census.setCache(new AtomCache(pdbDir, false));
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

	public static void main(String[] args) {
		final String pdbDir = args[0];
		final File censusFile = new File(args[1]);
		final int identityCutoff = Integer.parseInt(args[2]);
		int[] sunIds = new int[args.length - 2];
		for (int i = 3; i < args.length; i++) {
			sunIds[i - 3] = Integer.parseInt(args[i]);
		}
		buildDefault(pdbDir, censusFile, identityCutoff, sunIds);
	}

	public AstralScopDescriptionCensus(int maxThreads, int identityCutoff, int[] sunIds) {
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
