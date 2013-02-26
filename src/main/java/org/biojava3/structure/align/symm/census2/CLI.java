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
 * Created on 2013-02-20
 *
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
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopInstallation;

/**
 * Decent command-line interface to Census. Combines aspects of many of the other census subclasses.
 * @author dmyerstu
 */
public class CLI {

	static final Logger logger = Logger.getLogger(CensusJob.class.getPackage().getName());

	private static final String NEWLINE = "\n";

	private static Options getOptions() {
		Options options = new Options();
		options.addOption(OptionBuilder.hasArg(true).withDescription("The PDB directory. Defaults to the AtomCache default.").isRequired(false).create("pdb"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("The output file. Defaults to ./census.xml.").isRequired(false).create("file"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("The number of threads to use. Defaults to the number available minus 1.").isRequired(false).create("threads"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Ignore any existing work and start from scratch.").isRequired(false).create("restart"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Don't use Berkeley SCOP.").isRequired(false).create("oldscop"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Prefetch all PDB files.").isRequired(false).create("prefetch"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Write to file every n jobs.").isRequired(false).create("every"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Apply sequence clustering at the specified sequence identity, in decimal. A value of 0 (default) means no sequence clustering.").isRequired(false).create("clustering"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Run on only the specified space-seperated list of SCOP sun ids. Use either this or -superfamilies.").isRequired(false).create("sunids"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Run on only the specified space-seperated list of superfamily names. Use either this or -sunids. Defaults to SCOP classes A-F, \"46456 48724 51349 53931 56572 56835\"").isRequired(false).create("superfamilies"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Use only the specified number of entities from the sun ids or superfamilies selected. A value of -1 (default) means all.").isRequired(false).create("number"));
		options.addOption(OptionBuilder.hasArg(false).withDescription("Randomize the order in which the entities are chosen (from sun ids or superfamilies).").isRequired(false).create("randomize"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("The verbosity of debug statements printed to stdout.").isRequired(false).create("verbosity"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Run on only the nth fold. Special option for running on OSG.").isRequired(false).create("foldindex"));
		options.addOption(OptionBuilder.hasArg(true).withDescription("Run on only the nth superfamily. Special option for running on OSG.").isRequired(false).create("sfindex"));
		return options;
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
			
			ScopFactory.setScopDatabase(new BerkeleyScopInstallation());

			final String pdbDir = cmd.getOptionValue("pdb");
			final String censusFile = cmd.getOptionValue("file");
			final Integer nThreads = cmd.getOptionValue("threads")==null? null : Integer.parseInt(cmd.getOptionValue("threads"));
			final Integer writeEvery = cmd.getOptionValue("every")==null? null : Integer.parseInt(cmd.getOptionValue("every"));
			final Integer number = cmd.getOptionValue("number")==null? null : Integer.parseInt(cmd.getOptionValue("number"));
			final Integer clustering = cmd.getOptionValue("clustering")==null? null : Integer.parseInt(cmd.getOptionValue("clustering"));
			
			int[] sunIds = null;
			if (cmd.getOptionValue("sfindex") != null) {
				int index = Integer.parseInt(cmd.getOptionValue("sfindex"));
				sunIds = new int[] {sfByIndex(index)};
			} else if (cmd.getOptionValue("foldindex") != null) {
					int index = Integer.parseInt(cmd.getOptionValue("foldindex"));
					sunIds = new int[] {foldByIndex(index)};
			} else if (cmd.getOptionValue("sunids") != null) {
				String parts[] = cmd.getOptionValue("sunids").split(" ");
				sunIds = new int[parts.length];
				for (int i = 0; i < parts.length; i++) {
					sunIds[i] = Integer.parseInt(parts[i]);
				}
			}
			
			if (cmd.getOptionValue("verbosity") != null) Logger.getRootLogger().setLevel(Level.toLevel(cmd.getOptionValue("verbosity")));
			final String[] superfamilies = cmd.getOptionValue("superfamilies")==null? null : cmd.getOptionValue("superfamilies").split(" ");
			final boolean randomize = cmd.hasOption("randomize");
			final boolean restart = cmd.hasOption("restart");
			final boolean prefetch = cmd.hasOption("prefetch");
			final boolean oldScop = cmd.hasOption("oldscop");

			run(pdbDir, censusFile, nThreads, writeEvery, number, clustering, sunIds, superfamilies, randomize, restart, prefetch, oldScop);

		} catch (RuntimeException e) {
			printError(e);
		}
	}

	public static int foldByIndex(int index) {
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> allFolds = scop.getByCategory(ScopCategory.Fold);
		return allFolds.get(index).getSunID();
	}

	public static int sfByIndex(int index) {
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> allSfs = scop.getByCategory(ScopCategory.Superfamily);
		return allSfs.get(index).getSunID();
	}

	public static void run(final String pdbDir, final String censusFile, final Integer pNThreads, final Integer writeEvery, final Integer number,
			final Integer clustering, final int[] pSunIds, final String[] superfamilies, final boolean randomize, final boolean restart,
			boolean prefetch, final boolean oldScop) {

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
			protected Results getStartingResults() {
				if (restart) return new Results();
				return super.getStartingResults();
			}
			@Override
			protected List<ScopDomain> getDomains() {

				final ScopDatabase scop = ScopFactory.getSCOP();
				List<ScopDomain> domains = new ArrayList<ScopDomain>();

				// first we make a list of sun IDs to use
				List<Integer> sunIds = new ArrayList<Integer>();
				if (superfamilies == null) {
					int[] ppSunIds;
					if (pSunIds == null) {
						ppSunIds = new int[] {46456, 48724, 51349, 53931, 56572, 56835};
					} else {
						ppSunIds = pSunIds;
					}
					for (Integer sunId : ppSunIds) {
						sunIds.add(sunId);
					}
				} else {
					// get sun IDs from superfamily names
					List<ScopDescription> allSfs = scop.getByCategory(ScopCategory.Superfamily);
					for (ScopDescription superfamily : allSfs) {
						for (int i = 0; i < superfamilies.length; i++) {
							if (superfamilies[i].equals(superfamily.getClassificationId())) {
								sunIds.add(superfamily.getSunID());
								break;
							}
						}
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
					logger.info("Using sequence clustering at " + clustering + "%");
				}

				// okay, now create the final set
				for (int sunId : sunIds) {
					List<ScopDomain> putative = new ArrayList<ScopDomain>();
					ScopDescriptionCensus.getDomainsUnder(sunId, putative);
					if (randomize) { // randomize if we need to
						Collections.shuffle(putative);
						logger.debug("Taking " + (number==null? "all" : number) + " ids in random order");
					} else {
						logger.debug("Taking " + (number==null? "all" : number) + " ids in sequential order");
					}
					int numForId = 0; // keep track of the number we've added
					for (ScopDomain domain : putative) {
						if (clustering == null || clusterRepresentatives.contains(domain.getScopId())) {
							if (number != null && numForId >= number) break; // number==null means no limit
							domains.add(domain);
							numForId++;
						}
					}
				}
				logger.info("Found " + domains.size() + " domains");
				return domains;

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

	/**
	 * Prints an error message for {@code e} that shows causes and suppressed messages recursively. Just a little more
	 * useful than {@code e.printStackTrace()}.
	 * 
	 * @param e
	 */
	public static void printError(Exception e) {
		System.err.println(printError(e, ""));
	}

	private static void printUsage(String note, Options options) {
		if (note != null) System.out.println(note);
		HelpFormatter hf = new HelpFormatter();
		hf.printHelp("java -jar CensusCLI.jar", options);
	}

}
