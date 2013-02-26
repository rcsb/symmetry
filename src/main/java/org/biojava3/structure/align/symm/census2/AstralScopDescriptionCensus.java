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
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * A census that takes a list of sun Ids and a sequence clustering.
 * @author dmyerstu
 */
public class AstralScopDescriptionCensus extends ScopDescriptionCensus {

	private static final String SERVER_LOCATION = FarmJobParameters.DEFAULT_SERVER_URL;
	
	static final Logger logger = Logger.getLogger(CensusJob.class.getPackage().getName());
	private int[] sunIds;
	private int identityCutoff;

	static {
		BasicConfigurator.configure();
	}

	public static void buildDefault(String pdbDir, File censusFile, int identityCutoff, int[] sunIds) {
		try {
			ScopFactory.setScopDatabase(new BerkeleyScopInstallation());
			System.setProperty(AbstractUserArgumentProcessor.PDB_DIR, pdbDir); // okay, this doesn't appear to be
																				// working
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			if (maxThreads < 1) maxThreads = 1;
			AstralScopDescriptionCensus census = new AstralScopDescriptionCensus(maxThreads, identityCutoff, sunIds);
			census.setOutputWriter(censusFile);
			census.setCache(new AtomCache(pdbDir, false));
			census.run();
			System.out.println(census);
		} catch (RuntimeException e) {
			logger.fatal(e.getMessage(), e);
		}
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
				if (clusterRepresentatives.contains(domain.getScopId())) {
					domains.add(domain);
				}
			}
		}
		return domains;
	}

	public static Set<String> getClusterRepresentatives(int cutoff) {
		Set<String> names = JFatCatClient.getRepresentatives(SERVER_LOCATION, cutoff);
		if (names == null) throw new MissingResourceException("Could not retrieve representatives from " + SERVER_LOCATION, "JFatCatClient",  "representative chains");
		return names;
	}

}
