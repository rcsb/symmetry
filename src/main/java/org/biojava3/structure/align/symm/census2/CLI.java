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
package org.biojava3.structure.align.symm.census2;

import java.io.File;
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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopInstallation;
import org.biojava3.structure.align.symm.census2.AstralScopDescriptionCensus.AstralSet;

/**
 * Decent command-line interface to Census. Combines aspects of many of the other census subclasses.
 * 
 * @author dmyerstu
 */
public class CLI {

	private static final Logger logger = LogManager.getLogger(Census.class.getPackage().getName());

	private static final String NEWLINE = "\n";

	public static int foldByIndex(int index) {
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> allFolds = scop.getByCategory(ScopCategory.Fold);
		return allFolds.get(index).getSunID();
	}

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

			ScopDatabase scop = ScopFactory.getSCOP();
			if (!scop.getClass().getName().equals(BerkeleyScopInstallation.class.getName())) { // for efficiency
				ScopFactory.setScopDatabase(new BerkeleyScopInstallation());
			}

			final String pdbDir = cmd.getOptionValue("pdb");
			final String censusFile = cmd.getOptionValue("file");
			final Integer nThreads = cmd.getOptionValue("threads") == null ? null : Integer.parseInt(cmd
					.getOptionValue("threads"));
			final Integer writeEvery = cmd.getOptionValue("every") == null ? null : Integer.parseInt(cmd
					.getOptionValue("every"));
			final Integer number = cmd.getOptionValue("number") == null ? null : Integer.parseInt(cmd
					.getOptionValue("number"));
			// final Integer clustering = cmd.getOptionValue("clustering") == null ? null :
			// Integer.parseInt(cmd.getOptionValue("clustering"));
			final AstralSet clustering = cmd.getOptionValue("clustering") == null ? null : AstralSet.parse(cmd
					.getOptionValue("clustering"));

			int[] sunIds = null;
			if (cmd.getOptionValue("sfindex") != null) {
				if (cmd.hasOption("oldscop")) {
					printUsage(null, options);
					System.exit(-1);
				}
				int index = Integer.parseInt(cmd.getOptionValue("sfindex"));
				sunIds = new int[] { sfByIndex(index) };
			} else if (cmd.getOptionValue("foldindex") != null) {
				if (cmd.hasOption("oldscop")) {
					printUsage(null, options);
					System.exit(-1);
				}
				int index = Integer.parseInt(cmd.getOptionValue("foldindex"));
				sunIds = new int[] { foldByIndex(index) };
			} else if (cmd.getOptionValue("sunids") != null) {
				String parts[] = cmd.getOptionValue("sunids").split(" ");
				sunIds = new int[parts.length];
				for (int i = 0; i < parts.length; i++) {
					sunIds[i] = Integer.parseInt(parts[i]);
				}
			}

			final String[] superfamilies = cmd.getOptionValue("superfamilies") == null ? null : cmd.getOptionValue(
					"superfamilies").split(" ");
			final String[] folds = cmd.getOptionValue("folds") == null ? null : cmd.getOptionValue("folds").split(" ");
			final boolean randomize = cmd.hasOption("randomize");
			final boolean restart = cmd.hasOption("restart");
			final boolean prefetch = cmd.hasOption("prefetch");
			final boolean oldScop = cmd.hasOption("oldscop");

			final String sigClass = cmd.getOptionValue("sigclass");
			final String sigMethod = cmd.getOptionValue("sigmethod");

			run(pdbDir, censusFile, nThreads, writeEvery, number, clustering, sunIds, superfamilies, folds, randomize,
					restart, prefetch, oldScop, sigClass, sigMethod);

		} catch (RuntimeException e) {
			printError(e);
		}
	}

	/**
	 * Prints an error message for {@code e} that shows causes and suppressed messages recursively. Just a little more
	 * useful than {@code e.printStackTrace()}.
	 * 
	 * @param e
	 */
	public static void printError(Exception e) {
		System.err.println(printError(e, ""));
	}

	public static void run(final String pdbDir, final String censusFile, final Integer pNThreads,
			final Integer writeEvery, final Integer number, final AstralSet clustering, final int[] pSunIds,
			final String[] superfamilies, final String[] folds, final boolean randomize, final boolean restart,
			boolean prefetch, final boolean oldScop, String sigClass, String sigMethod) {

		// get a significance object
		final Significance sig;
		if (sigClass != null) {
			if (sigMethod != null) {
				sig = SignificanceFactory.fromMethod(sigClass, sigMethod);
			} else {
				sig = SignificanceFactory.fromClass(sigClass);
			}
		} else {
			if (sigMethod != null) {
				sig = SignificanceFactory.fromMethod(SignificanceFactory.class.getName(), sigMethod);
			} else {
				sig = SignificanceFactory.forCensus();
			}
		}

		// set the number of threads
		final int nThreads;
		if (pNThreads == null) {
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			if (maxThreads < 1) maxThreads = 1;
			nThreads = maxThreads;
		} else {
			nThreads = pNThreads;
		}
		logger.info("Using " + nThreads + " threads");

		Census census = new Census(nThreads) {
			@Override
			protected List<ScopDomain> getDomains() {

				// first we make a list of sun IDs to use
				List<Integer> sunIds = new ArrayList<Integer>();
				if (pSunIds != null) {
					for (Integer pSunId : pSunIds) {
						sunIds.add(pSunId);
					}
				}
				sunIds.addAll(getSunIds(superfamilies, ScopCategory.Superfamily));
				sunIds.addAll(getSunIds(folds, ScopCategory.Fold));
				if (sunIds.isEmpty()) {
					int[] ppSunIds = new int[] { 46456, 48724, 51349, 53931, 56572, 56835 };
					for (Integer sunId : ppSunIds) {
						sunIds.add(sunId);
					}
				}

				// print the sun IDs we're using
				StringBuilder sb = new StringBuilder();
				for (int i = 0; i < sunIds.size(); i++) {
					sb.append(sunIds.get(i));
					if (i < sunIds.size() - 1) sb.append(", ");
				}
				logger.info("Using sun IDs " + sb.toString());

				// get sequence clusters
				Set<String> clusterRepresentatives = null;
				if (clustering != null) {
					clusterRepresentatives = AstralScopDescriptionCensus.getClusterRepresentatives(clustering);
					logger.info("Using sequence clustering at " + clustering.getId());
					logger.info("Number of clusters: " + clusterRepresentatives.size());
				}

				// okay, now create the final set
				List<ScopDomain> domains = new ArrayList<ScopDomain>();
				for (int sunId : sunIds) {

					// first, let putative contain all the domains under our sun id
					List<ScopDomain> putative = new ArrayList<ScopDomain>();
					ScopDescriptionCensus.getDomainsUnder(sunId, putative);

					// randomize if we need to
					if (randomize) {
						Collections.shuffle(putative);
						logger.debug("Taking " + (number == null ? "all" : number) + " ids in random order out of "
								+ putative.size());
					} else {
						logger.debug("Taking " + (number == null ? "all" : number) + " ids in sequential order out of "
								+ putative.size());
					}

					// keep track of the number we've added
					// we don't want to add more than the number we're allowed to
					int numForId = 0;

					for (ScopDomain domain : putative) {

						boolean contains = clustering == null;
						if (!contains) {
							contains = clusterRepresentatives.contains(domain.getScopId());
						}
						// if (!contains) {
						// contains = PdbClusteringScopDescriptionCensus.isDomainOverChain(domain,
						// clusterRepresentatives);
						// }
						// logger.debug("Contains " + domain.getScopId());

						if (contains) { // if we want to include the domain
							// don't add more than we're allowed to
							if (number != null && numForId >= number) break; // number == null means no limit
							domains.add(domain);
							numForId++;
						}
					}
				}
				logger.info("Found " + domains.size() + " domains");
				return domains;

			}

			@Override
			protected Significance getSignificance() {
				return sig;
			}

			@Override
			protected Results getStartingResults() {
				if (restart) return new Results();
				return super.getStartingResults();
			}
		};

		// set final options
		if (!oldScop) {
			ScopFactory.setScopDatabase(new BerkeleyScopInstallation());
		} else {
			ScopFactory.setScopDatabase(new ScopInstallation());
		}
		if (writeEvery != null) census.setPrintFrequency(writeEvery);
		census.setDoPrefetch(prefetch);
		if (censusFile != null) {
			census.setOutputWriter(new File(censusFile));
		} else {
			census.setOutputWriter(new File("census.xml"));
		}
		if (pdbDir == null) {
			census.setCache(new AtomCache());
		} else {
			census.setCache(new AtomCache(pdbDir, false));
		}

		// now run
		census.run();
		System.out.println(census);
	}

	public static int sfByIndex(int index) {
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> allSfs = scop.getByCategory(ScopCategory.Superfamily);
		return allSfs.get(index).getSunID();
	}

	private static Options getOptions() {
		Options options = new Options();
		options.addOption(OptionBuilder.hasArg(true)
				.withDescription("The PDB directory. Defaults to the AtomCache default.").isRequired(false)
				.create("pdb"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("The output file. Defaults to ./census.xml.")
				.isRequired(false).create("file"));
		options.addOption(OptionBuilder.hasArg(true)
				.withDescription("The number of threads to use. Defaults to the number available minus 1.")
				.isRequired(false).create("threads"));
		options.addOption(OptionBuilder.hasArg(false)
				.withDescription("Ignore any existing work and start from scratch.").isRequired(false)
				.create("restart"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Don't use Berkeley SCOP.").isRequired(false)
				.create("oldscop"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Prefetch all PDB files.").isRequired(false)
				.create("prefetch"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Write to file every n jobs.").isRequired(false)
				.create("every"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Apply sequence clustering at the specified sequence identity, in decimal. A value of 0 (default) means no sequence clustering.")
				.isRequired(false).create("clustering"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Run on only the specified space-seperated list of SCOP sun ids. Defaults to SCOP classes A-F, \"46456 48724 51349 53931 56572 56835\"")
				.isRequired(false).create("sunids"));
		options.addOption(OptionBuilder.hasArg(true)
				.withDescription("Run on only the specified space-seperated list of superfamily names.")
				.isRequired(false).create("superfamilies"));
		options.addOption(OptionBuilder.hasArg(true)
				.withDescription("Run on only the specified space-seperated list of fold names.").isRequired(false)
				.create("folds"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Use only the specified number of entities from the sun ids or superfamilies selected. A value of -1 (default) means all. This option does not attempt to select diverse domains from the sets (will fix). Use with -randomize.")
				.isRequired(false).create("number"));
		options.addOption(OptionBuilder
				.hasArg(false)
				.withDescription(
						"Randomize the order in which the entities are chosen (from sun ids or superfamilies).")
				.isRequired(false).create("randomize"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Run on only the nth fold. Special option for running on OSG. Does not work with -oldscop.")
				.isRequired(false).create("foldindex"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"Run on only the nth superfamily. Special option for running on OSG. Does not work with -oldscop.")
				.isRequired(false).create("sfindex"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"A fully-qualified class name for a Significance object. If sigmethod is also set, calls that factory method; otherwise, calls a constructor.")
				.isRequired(false).create("sigclass"));
		options.addOption(OptionBuilder
				.hasArg(true)
				.withDescription(
						"The name of a factory method that returns a Significance object. If sigclass is also set, expects the factory method to be in that class; otherwise, checks in SignificanceFactory.")
				.isRequired(false).create("sigmethod"));
		return options;
	}

	private static List<Integer> getSunIds(String[] superfamilies, ScopCategory category) {
		List<Integer> sunIds = new ArrayList<Integer>();
		if (superfamilies == null) return sunIds;
		final ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> allSfs = scop.getByCategory(category);
		for (ScopDescription superfamily : allSfs) {
			for (String superfamilie : superfamilies) {
				if (superfamilie.equals(superfamily.getClassificationId())) {
					sunIds.add(superfamily.getSunID());
					break;
				}
			}
		}
		return sunIds;
	}

	/**
	 * @see #printError(Exception)
	 */
	private static String printError(Exception e, String tabs) {
		StringBuilder sb = new StringBuilder();
		Throwable prime = e;
		while (prime != null) {
			if (tabs.length() > 0) sb.append(tabs + "Cause:" + NEWLINE);
			sb.append(tabs + prime.getClass().getSimpleName() + NEWLINE);
			if (prime.getMessage() != null) sb.append(tabs + prime.getMessage() + NEWLINE);
			if (prime instanceof Exception) {
				StackTraceElement[] trace = ((Exception) prime).getStackTrace();
				for (StackTraceElement element : trace) {
					sb.append(tabs + element.toString() + NEWLINE);
				}
			}
			prime = prime.getCause();
			tabs += "\t";
			sb.append(NEWLINE);
		}
		return sb.toString();
	}

	private static void printUsage(String note, Options options) {
		if (note != null) System.out.println(note);
		HelpFormatter hf = new HelpFormatter();
		hf.printHelp("java -jar CensusCLI.jar", options);
	}

}
