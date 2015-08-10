/*
 * BioJava development code
 * 
 * This code may be freely distributed and modified under the terms of the GNU Lesser General Public Licence. This
 * should be distributed with the code. If you do not have a copy, see:
 * 
 * http://www.gnu.org/copyleft/lesser.html
 * 
 * Copyright for this code is held jointly by the individual authors. These should be listed in @author doc comments.
 * 
 * For more information on the BioJava project and its aims, or to join the biojava-l mailing list, visit the home page
 * at:
 * 
 * http://www.biojava.org/
 * 
 * Created on 2013-02-20
 */
package org.biojava.nbio.structure.align.symm.census3.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.scop.Astral;
import org.biojava.nbio.structure.scop.Astral.AstralSet;
import org.biojava.nbio.structure.scop.ScopCategory;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.census3.CensusResultList;
import org.biojava.nbio.structure.align.symm.census3.representatives.ScopSupport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Decent command-line interface to Census. Combines aspects of many of the other census subclasses.
 * 
 * @author dmyersturnbull
 */
public class CensusCLI {

	private final static Logger logger = LoggerFactory.getLogger(CensusCLI.class);

	/**
	 * See {@link #run(String, String, String, String, String, PrintStream)}.
	 * 
	 * @param args
	 * @throws ServiceException
	 * @throws IOException
	 * @throws NameException
	 */
	public static void main(String[] args) {

		try {

			Options options = getOptions();
			CommandLineParser parser = new GnuParser();
			CommandLine cmd;
			try {
				cmd = parser.parse(options, args);
			} catch (ParseException e) {
				printUsage(null, options);
				return;
			}

			// set SCOP version
			String scopVersion = cmd.getOptionValue("scopversion");
			if (scopVersion != null && !scopVersion.isEmpty())
				ScopFactory.setScopDatabase(scopVersion);

			final String pdbDir = cmd.getOptionValue("pdbdir");
			final String censusFile = cmd.getOptionValue("file");
			final Integer nThreads = cmd.getOptionValue("threads") == null ? null : Integer.parseInt(cmd
					.getOptionValue("threads"));
			final Integer writeEvery = cmd.getOptionValue("every") == null ? null : Integer.parseInt(cmd
					.getOptionValue("every"));
			final Integer number = cmd.getOptionValue("maxreps") == null ? null : Integer.parseInt(cmd
					.getOptionValue("maxreps"));
			final AstralSet clustering = cmd.getOptionValue("clustering") == null ? null : AstralSet.parse(cmd
					.getOptionValue("clustering"));

			int[] sunIds = null;
			if (cmd.getOptionValue("sfindex") != null) {
				int index = Integer.parseInt(cmd.getOptionValue("sfindex"));
				sunIds = new int[] { ScopSupport.getInstance().getByIndex(ScopCategory.Superfamily, index).getSunID() };
			} else if (cmd.getOptionValue("foldindex") != null) {
				int index = Integer.parseInt(cmd.getOptionValue("foldindex"));
				sunIds = new int[] { ScopSupport.getInstance().getByIndex(ScopCategory.Fold, index).getSunID() };
			} else if (cmd.getOptionValue("sunids") != null) {
				String parts[] = cmd.getOptionValue("sunids").split(" ");
				sunIds = new int[parts.length];
				for (int i = 0; i < parts.length; i++) {
					sunIds[i] = Integer.parseInt(parts[i]);
				}
			} else if (cmd.getOptionValue("names") != null) {
				// get sun ids from the SCOP ids
				// this is pretty stupid, since we'll do the opposite in a bit
				ScopDatabase scop = ScopFactory.getSCOP();
				BufferedReader br = null;
				try {
					br = new BufferedReader(new FileReader(new File(cmd.getOptionValue("names"))));
					String line = "";
					List<Integer> list = new ArrayList<Integer>();
					while ((line = br.readLine()) != null) {
						ScopDomain domain = scop.getDomainByScopID(line.trim());
						list.add(domain.getSunid());
					}
					sunIds = new int[list.size()];
					for (int i = 0; i < list.size(); i++) {
						sunIds[i] = list.get(i);
					}
				} catch (IOException e) {
					throw new RuntimeException("Couldn't parse list of SCOP Ids", e);
				} finally {
					if (br != null) {
						try {
							br.close();
						} catch (IOException e) {
							logger.warn("Couldn't close file " + cmd.getOptionValue("names"), e);
						}
					}
				}
			}
			final int start = cmd.getOptionValue("start") == null ? 0 : Integer.parseInt(cmd
					.getOptionValue("start"));
			final int stop = cmd.getOptionValue("stop") == null ? Integer.MAX_VALUE : Integer.parseInt(cmd
					.getOptionValue("stop"));

			final String[] classifications = cmd.getOptionValue("classids") == null ? null : cmd.getOptionValue(
					"classids").split(" ");
			final boolean randomize = cmd.hasOption("randomize");
			final boolean restart = cmd.hasOption("restart");
			final boolean prefetch = cmd.hasOption("prefetch");
			final boolean noMapping = cmd.hasOption("nomapping");
			final boolean diverse = cmd.hasOption("diverse");
			final boolean allProteins = cmd.hasOption("allproteins");
			final boolean reverse = cmd.hasOption("reverse");

			run(pdbDir, censusFile, nThreads, writeEvery, number, clustering, sunIds, classifications, randomize,
					restart, prefetch, noMapping, scopVersion, diverse, allProteins, start, stop, reverse);

		} catch (RuntimeException e) {
			logger.error(e.getMessage(),e);
		}
	}

	/**
	 * Returns the best number of threads.
	 */
	private static int getNThreads(Integer inputNThreads) {
		// set the number of threads
		final int nThreads;
		if (inputNThreads == null) {
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			if (maxThreads < 1) maxThreads = 1;
			nThreads = maxThreads;
		} else {
			nThreads = inputNThreads;
		}
		return nThreads;
	}

	/**
	 * Gets an <em>initial</em> list of sun ids to use, without clustering.
	 */
	private static List<Integer> getSunIdsToUse(int[] inputSunIds, String[] inputClassIds, AstralSet clustering) {
		List<Integer> sunIds = new ArrayList<Integer>();
		if (inputSunIds != null) {
			for (Integer sunId : inputSunIds) {
				sunIds.add(sunId);
			}
		}
		if (sunIds.isEmpty()) {
			for (Integer sunId : ScopSupport.TRUE_SCOP_CLASSES) {
				sunIds.add(sunId);
			}
		}
		return sunIds;
	}

	/**
	 * Returns a <em>clustered and truncated</em> list of ScopDomains from the list of sun ids.
	 */
	private static List<ScopDomain> getDomainsFromSunIds(List<Integer> sunIds, AstralSet clustering, Integer limit, boolean diverse, boolean allProteins, boolean randomize) {

		// get sequence clusters
		Set<String> clusterRepresentatives = null;
		if (clustering != null) {
			clusterRepresentatives = Astral.getRepresentatives(clustering);
			logger.info("Using sequence clustering at " + clustering.getId());
			logger.info("Number of clusters: " + clusterRepresentatives.size());
		}

		// okay, now create the final set
		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		for (int sunId : sunIds) {

			// first, let putative contain all the domains under our sun id
			List<ScopDomain> putative = new ArrayList<ScopDomain>();
			if (diverse) {
				ScopSupport.getInstance().getDomainsUnder(sunId, putative, limit, allProteins);
			} else {
				ScopSupport.getInstance().getAllDomainsUnder(sunId, putative, allProteins);
			}

			// randomize if we need to
			if (randomize) {
				Collections.shuffle(putative);
				logger.debug("Taking " + (limit == null ? "all" : limit) + " ids in random order out of "
						+ putative.size());
			} else {
				logger.debug("Taking " + (limit == null ? "all" : limit) + " ids in sequential order out of "
						+ putative.size());
			}

			// now handle clustering and limited number
			// keep track of the number we've added
			// we don't want to add more than the number we're allowed to
			int numForId = 0;

			for (ScopDomain domain : putative) {

				boolean contains = clustering == null;
				if (!contains) {
					contains = clusterRepresentatives.contains(domain.getScopId());
				}
				if (contains) { // if we want to include the domain
					// don't add more than we're allowed to
					if (limit != null && numForId >= limit) break; // number == null means no limit
					domains.add(domain);
					numForId++;
				}
			}
		}

		return domains;

	}

	public static void run(final String pdbDir, final String censusFile, final Integer inputNThreads,
			final Integer writeEvery, final Integer number, final AstralSet clustering, final int[] inputSunIds,
			final String[] inputClassIds, final boolean randomize, final boolean restart,
			boolean prefetch, final boolean noMapping, String scopVersion, final boolean diverse, final boolean allProteins, final int start, final int stop, final boolean reverse) {

		Census census;

		{

			final int nThreads = getNThreads(inputNThreads);

			logger.info("Using " + nThreads + " threads");

			census = new Census(nThreads) {
				@Override
				protected List<ScopDomain> getDomains() {

					// get list of sun Ids
					List<Integer> sunIds = getSunIdsToUse(inputSunIds, inputClassIds, clustering);

					// print the sun IDs we're using
					StringBuilder sb = new StringBuilder();
					for (int i = 0; i < sunIds.size(); i++) {
						sb.append(sunIds.get(i));
						if (i < sunIds.size() - 1) sb.append(", ");
					}
					logger.info("Using sun IDs " + sb.toString());

					List<ScopDomain> domains = getDomainsFromSunIds(sunIds, clustering, number, diverse, allProteins, randomize);
					List<ScopDomain> restrictedList = new ArrayList<ScopDomain>(Math.min(domains.size(), stop - start));
					for (int i = start; i < Math.min(domains.size(), stop); i++) {
						restrictedList.add(domains.get(i));
					}
					logger.info("Found " + domains.size() + " domains and using " + restrictedList.size());
					if (reverse) {
						Collections.reverse(restrictedList);
					}
					return restrictedList;

				}

				@Override
				protected CensusResultList getStartingResults() {
					if (restart) return new CensusResultList();
					return super.getStartingResults();
				}
			};

			// set PDB dir
			// this actually gets called first
			if (pdbDir == null) {
				AtomCache cache = new AtomCache();
				cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
				census.setCache(cache);
			} else {
				AtomCache cache = new AtomCache(pdbDir,pdbDir);
				cache.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
				cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
				census.setCache(cache);
				System.setProperty(UserConfiguration.PDB_DIR, pdbDir);
			}

			// set final options
			if (writeEvery != null) census.setPrintFrequency(writeEvery);
			census.setDoPrefetch(prefetch);
			if (censusFile != null) {
				census.setOutputWriter(new File(censusFile));
			} else {
				census.setOutputWriter(new File("census.xml"));
			}
			census.setRecordAlignmentMapping(!noMapping);

		}

		// now run
		ScopSupport.nullifyInstance();
		census.run();
		System.out.println(census);
	}

	@SuppressWarnings("static-access")
	private static Options getOptions() {
		Options options = new Options();
		options.addOption(OptionBuilder.hasArg(true)
				.withDescription("The PDB directory. Defaults to the Biojava AtomCache default, which is typically the operating system's temporary directory.").isRequired(false)
				.create("pdbdir"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("The output file. Defaults to \"./census.xml\".")
				.isRequired(false).create("file"));
		options.addOption(OptionBuilder.hasArg(true)
				.withDescription("The number of threads to use. Defaults to the number available minus 1.")
				.isRequired(false).create("threads"));
		options.addOption(OptionBuilder.hasArg(false)
				.withDescription("Ignore any existing work and start from scratch. Otherwise, the census will read the output file and not repeat any jobs.").isRequired(false)
				.create("restart"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Use the specified SCOP version; defaults to the most recent.").isRequired(false)
				.create("scopversion"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Prefetch all PDB files before starting.").isRequired(false)
				.create("prefetch"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Do not record the alignment mapping in the XML. The mapping can be used to reconstruct an AFPChain quickly.").isRequired(false)
				.create("nomapping"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Use all \"px\"s under every domain. If not set, uses only the first \"px\" of each domain.").isRequired(false)
				.create("allproteins"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Reverses the direction in which the jobs are run.").isRequired(false)
				.create("reverse"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("If -maxreps is less than the total number of domains for a SCOP category, tries to spread the selection over the category. Currently only works with SCOP superfamilies.").isRequired(false)
				.create("diverse"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Write to the file every n jobs.").isRequired(false)
				.create("every"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Start writing on the nth domain selected. Defaults to 0.").isRequired(false)
				.create("start"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Stop writing after the nth domains selected. Defaults to +infinity.").isRequired(false)
				.create("stop"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Use a file containing a line-by-line list of SCOP Ids to run on.").isRequired(false)
				.create("names"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Apply sequence clustering using ASTRAL at the specified sequence identity, in decimal. A value of 0 (default) means no sequence clustering.")
						.isRequired(false).create("clustering"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Run on only the specified space-seperated list of SCOP sun ids. Defaults to SCOP classes A-F, \"46456 48724 51349 53931 56572 56835\".")
						.isRequired(false).create("sunids"));
		options.addOption(OptionBuilder.hasArg(true)
				.withDescription("Run on only the specified space-seperated list of classification identifiers, in addition to -sunids.")
				.isRequired(false).create("classids"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Use at most the specified number of entities from the sun ids or superfamilies selected. A value of -1 (default) means all. This option does not attempt to select diverse domains from the sets (will fix). Use with -randomize.")
						.isRequired(false).create("maxreps"));
		options.addOption(OptionBuilder
				.hasArg(false)
				.withDescription(
						"Randomize the order in which the entities are chosen (from sun ids or superfamilies).")
						.isRequired(false).create("randomize"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Run on only the \"nth\" fold, where the ordering is arbitrary. Historically used for running on Open Science Grid.")
						.isRequired(false).create("foldindex"));
		return options;
	}

	private static void printUsage(String note, Options options) {
		if (note != null) System.out.println(note);
		HelpFormatter hf = new HelpFormatter();
		hf.printHelp("java -jar CensusCLI.jar", options);
	}

}
