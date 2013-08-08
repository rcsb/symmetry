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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.representatives.ScopSupport;

/**
 * A census that takes a list of sun Ids and a sequence clustering.
 * 
 * @author dmyerstu
 */
public class AstralScopDescriptionCensus extends ScopDescriptionCensus {

	private static final Logger logger = LogManager.getLogger(Census.class.getPackage().getName());

	public static enum AstralSet {
		FORTY_175A("1.75A_40", "http://scop.berkeley.edu/downloads/scopseq-1.75A/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75A.fa"),
		NINETY_FIVE_175A("1.75A_95", "http://scop.berkeley.edu/downloads/scopseq-1.75A/astral-scopdom-seqres-gd-sel-gs-bib-90-1.75A.fa"),
		FORTY_175B("1.75B_95", "http://scop.berkeley.edu/downloads/scopseq-1.75B/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75B.fa"),
		NINETY_FIVE_175B("1.75B_95", "http://scop.berkeley.edu/downloads/scopseq-1.75B/astral-scopdom-seqres-gd-sel-gs-bib-95-1.75B.fa"),
		FORTY_175("1.75_95", "http://scop.berkeley.edu/downloads/scopseq-1.75/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"),
		NINETY_FIVE_175("1.75_95", "http://scop.berkeley.edu/downloads/scopseq-1.75/astral-scopdom-seqres-gd-sel-gs-bib-95-1.75.fa");
		private String url;
		private String id;
		public static AstralSet parse(String str) {
			for (AstralSet c : AstralSet.class.getEnumConstants()) {
				if (c.getId().equals(str)) return c;
			}
			throw new IllegalArgumentException("No ASTRAL set with id " + str);
		}
		AstralSet(String id, String url) {
			this.url = url;
			this.id = id;
		}
		String getId() {
			return id;
		}
		String getUrl() {
			return url;
		}
	}
	
	public static void buildDefault(File censusFile, AstralSet astral, int[] sunIds) {
		try {
			ScopFactory.setScopDatabase(ScopFactory.getSCOP(ScopFactory.VERSION_1_75A));
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			AstralScopDescriptionCensus census = new AstralScopDescriptionCensus(maxThreads, astral, sunIds);
			census.setOutputWriter(censusFile);
			census.setCache(new AtomCache());
			census.run();
			System.out.println(census);
		} catch (RuntimeException e) {
			logger.fatal(e.getMessage(), e);
		}
	}

	public static Set<String> getClusterRepresentatives(AstralSet cutoff) {
		
		URL url;
		try {
			url = new URL(cutoff.getUrl());
		} catch (MalformedURLException e) {
			throw new RuntimeException("The URL was invalid!", e);
		}
		
		Set<String> names = new TreeSet<String>();
		
		try {
			
			BufferedReader br = new BufferedReader(new InputStreamReader(url.openStream()));

			logger.info("Reading ASTRAL file...");
			
			String line = "";
			int i = 0;
			while ((line = br.readLine()) != null) {
				if (line.startsWith(">")) {
					try {
						String scopId = line.split("\\s")[0].substring(1);
						names.add(scopId);
						if (i % 1000 == 0) {
							logger.debug("Reading ASTRAL line for " + scopId);
						}
						i++;
					} catch (RuntimeException e) {
						logger.error("Couldn't read line " + line, e);
					}
				}
			}

			br.close();
			
		} catch (IOException e) {
			throw new RuntimeException("Couldn't read the input stream to " + url.getPath(), e);
		}

		return names;
	}

	public static void main(String[] args) {
		if (args.length < 3) {
			System.err.println("Usage: " + AstralScopDescriptionCensus.class.getSimpleName() + " output-census-file astral-version-id sun-id-1 [sun-id-2 sun-id-3 ...]");
			return;
		}
		final File censusFile = new File(args[0]);
		final AstralSet astral = AstralSet.parse(args[1]);
		int[] sunIds = new int[args.length - 2];
		StringBuilder sb = new StringBuilder();
		for (int i = 2; i < args.length; i++) {
			Integer sunId = ScopSupport.getInstance().getSunId(args[i]);
			if (sunId == null) throw new IllegalArgumentException("Couldn't find " + args[i]);
			sb.append(sunId + ",");
			sunIds[i - 2] = sunId;
		}
		logger.info("Running on " + sb.toString().substring(0, sb.toString().length()-1));
		buildDefault(censusFile, astral, sunIds);
	}

	private AstralSet astral;
	
	public AstralScopDescriptionCensus(int maxThreads, AstralSet astral, int[] sunIds) {
		super(maxThreads, sunIds);
		this.astral = astral;
	}

	@Override
	protected List<ScopDomain> getDomains() {

		final Set<String> clusterRepresentatives = getClusterRepresentatives(astral);
		List<ScopDomain> domains = new ArrayList<ScopDomain>();

		if (sunIds == null) throw new RuntimeException("WHAT?!");

		for (int sunId : sunIds) {
			
			List<ScopDomain> putative = new ArrayList<ScopDomain>();
			ScopDescriptionCensus.getDomainsUnder(sunId, putative);

			for (ScopDomain domain : putative) {
				if (clusterRepresentatives.contains(domain.getScopId())) domains.add(domain);
			}

		}

		return domains;
	}

}
