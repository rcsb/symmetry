package main;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.CliTools;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.io.StructureFiletype;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import workers.QuatSymmWorker;
import writers.QuatSymmFastaWriter;
import writers.QuatSymmStatsWriter;
import writers.QuatSymmWriter;

/**
 * Main executable for running the Quaternary Symmetry detection. Run with -h
 * for usage help.
 *
 * @author Aleix Lafita
 *
 */
public class QuatSymmMain {

	private static final Logger logger = LoggerFactory
			.getLogger(QuatSymmMain.class);

	public static void main(String[] args) throws InterruptedException {

		// Begin argument parsing
		final String usage = "[OPTIONS] [structures...]";
		final String header = "Determine the stoichiometry and symmetry of an assembly.";
		Options options = getOptions();
		CommandLineParser parser = new DefaultParser();
		HelpFormatter help = new HelpFormatter();
		help.setOptionComparator(null); // prevent option sorting

		final CommandLine cli;
		try {
			cli = parser.parse(options, args, false);
		} catch (ParseException e) {
			logger.error("Error: " + e.getMessage());
			help.printHelp(usage, header, options, "");
			System.exit(1);
			return;
		}

		args = cli.getArgs();

		// help
		if (cli.hasOption("help")) {
			help.printHelp(usage, header, options, "");
			System.exit(0);
			return;
		}

		// version
		if (cli.hasOption("version")) {
			String version = QuatSymmMain.class.getPackage()
					.getImplementationVersion();
			if (version == null || version.isEmpty()) {
				version = "(custom version)";
			}
			System.out.println("QuatSymm " + version);
			System.exit(0);
			return;
		}

		// input structures
		List<String> names = new ArrayList<String>();
		if (cli.hasOption("input")) {
			// read from file
			try {
				names = parseInputStructures(cli.getOptionValue("input"));
			} catch (FileNotFoundException e) {
				logger.error("Error: File not found: "
						+ cli.getOptionValue("input"));
				System.exit(1);
				return;
			}
			// append cli arguments
			names.addAll(Arrays.asList(args));
		} else {
			if (args.length == 0) {
				// No structures given; print help and return
				help.printHelp(usage, header, options, "");
				System.exit(0);
				return;
			} else {
				// take names from the command line arguments
				names = Arrays.asList(args);
			}
		}

		// Show jmol?
		// Default to false with --input or with >=10 structures
		boolean show3d = !cli.hasOption("input") && names.size() < 10;
		if (cli.hasOption("noshow3d")) {
			show3d = false;
		}
		if (cli.hasOption("show3d")) {
			show3d = true;
		}

		// AtomCache options
		String pdbFilePath = null;
		if (cli.hasOption("pdbfilepath")) {
			pdbFilePath = cli.getOptionValue("pdbfilepath");
			pdbFilePath = FileDownloadUtils.expandUserHome(pdbFilePath);
		}

		// Logger control
		if (cli.hasOption("verbose")) {
			// Note that this bypasses SLF4J abstractions
			LoggerContext ctx = (LoggerContext) LogManager.getContext(false);
			Configuration config = ctx.getConfiguration();
			LoggerConfig loggerConfig = config
					.getLoggerConfig(LogManager.ROOT_LOGGER_NAME);
			loggerConfig.setLevel(Level.INFO);
			ctx.updateLoggers(); // This causes all Loggers to refetch
									// information from their LoggerConfig.
		}

		// Output formats
		List<QuatSymmWriter> writers = new ArrayList<QuatSymmWriter>();

		if (cli.hasOption("stats")) {
			String filename = cli.getOptionValue("stats");
			if (filename == null || filename.isEmpty())
				filename = "-"; // standard out
			try {
				writers.add(new QuatSymmStatsWriter(filename));
			} catch (IOException e) {
				logger.error("Error: Ignoring file " + filename + ".");
				logger.error(e.getMessage());
			}
		}

		if (cli.hasOption("fasta")) {
			String filename = cli.getOptionValue("fasta");
			if(filename == null || filename.isEmpty())
				filename = "-"; // standard out
			try {
				writers.add(new QuatSymmFastaWriter(filename));
			} catch (IOException e) {
				logger.error("Error: Ignoring file " + filename + ".");
				logger.error(e.getMessage());
			}
		}

		// Default Writer
		if (writers.isEmpty() && !cli.hasOption("noverbose")) {
			try {
				writers.add(new QuatSymmStatsWriter("-"));
			} catch (IOException e) {
				logger.error(e.getMessage());
			}
		}

		// Multithreading
		Integer threads = Runtime.getRuntime().availableProcessors();
		if (cli.hasOption("threads")) {
			threads = new Integer(cli.getOptionValue("threads"));
			if (threads < 1) {
				threads = 1;
			}
		}

		// Subunit Clustering parameters
		SubunitClustererParameters cparams = new SubunitClustererParameters();

		if (cli.hasOption("minSeqLen")) {
			String value = cli.getOptionValue("minSeqLen");
			try {
				int minSeqLen = Integer.parseInt(value);
				if (minSeqLen < 1) {
					logger.error("Invalid minSeqLen: " + minSeqLen);
					System.exit(1);
				}
				cparams.setMinimumSequenceLength(minSeqLen);
			} catch (NumberFormatException e) {
				logger.error("Invalid minSeqLen: " + value);
				System.exit(1);
			}
		}
		if (cli.hasOption("minSeqId")) {
			String value = cli.getOptionValue("minSeqId");
			try {
				double minSeqId = Double.parseDouble(value);
				if (minSeqId < 0 || minSeqId > 1) {
					logger.error("Invalid minSeqId: " + minSeqId);
					System.exit(1);
				}
				cparams.setSequenceIdentityThreshold(minSeqId);
			} catch (NumberFormatException e) {
				logger.error("Invalid minSeqId: " + value);
				System.exit(1);
			}
		}
		if (cli.hasOption("minSequenceCoverage")) {
			String value = cli.getOptionValue("minSequenceCoverage");
			try {
				double minSequenceCoverage = Double.parseDouble(value);
				if (minSequenceCoverage < 0 || minSequenceCoverage > 1) {
					logger.error("Invalid minSequenceCoverage: " + minSequenceCoverage);
					System.exit(1);
				}
				cparams.setSequenceCoverageThreshold(minSequenceCoverage);
			} catch (NumberFormatException e) {
				logger.error("Invalid minSequenceCoverage: " + value);
				System.exit(1);
			}
		}
		if (cli.hasOption("minStructureCoverage")) {
			String value = cli.getOptionValue("minStructureCoverage");
			try {
				double minStructureCoverage = Double.parseDouble(value);
				if (minStructureCoverage < 0 || minStructureCoverage > 1) {
					logger.error("Invalid minStructureCoverage: " + minStructureCoverage);
					System.exit(1);
				}
				cparams.setStructureCoverageThreshold(minStructureCoverage);
			} catch (NumberFormatException e) {
				logger.error("Invalid minStructureCoverage: " + value);
				System.exit(1);
			}
		}
		if (cli.hasOption("maxClustRmsd")) {
			String value = cli.getOptionValue("maxClustRmsd");
			try {
				double maxClustRmsd = Double.parseDouble(value);
				if (maxClustRmsd < 0) {
					logger.error("Invalid maxClustRmsd: " + maxClustRmsd);
					System.exit(1);
				}
				cparams.setRMSDThreshold(maxClustRmsd);
			} catch (NumberFormatException e) {
				logger.error("Invalid maxClustRmsd: " + value);
				System.exit(1);
			}
		}
		if (cli.hasOption("clustMethod")) {
			String value = cli.getOptionValue("clustMethod");
			SubunitClustererMethod clustMethod;
			try {
				clustMethod = SubunitClustererMethod.valueOf(value
						.toUpperCase());
				cparams.setClustererMethod(clustMethod);
			} catch (IllegalArgumentException e) {
				// give up
				logger.error("Illegal subunit clusterer method. Requires on of "
						+ CliTools
								.getEnumValuesAsString(SubunitClustererMethod.class));
				System.exit(1);
			}
		}

		// Quaternary Symmetry Parameters
		QuatSymmetryParameters sparams = new QuatSymmetryParameters();

		if (cli.hasOption("maxSymmRmsd")) {
			String value = cli.getOptionValue("maxSymmRmsd");
			try {
				double maxSymmRmsd = Double.parseDouble(value);
				if (maxSymmRmsd < 0) {
					logger.error("Invalid maxSymmRmsd: " + maxSymmRmsd);
					System.exit(1);
				}
				sparams.setRmsdThreshold(maxSymmRmsd);
			} catch (NumberFormatException e) {
				logger.error("Invalid maxSymmRmsd: " + value);
				System.exit(1);
			}
		}

		// Done parsing arguments

		// Write the headers of the files
		for (QuatSymmWriter writer : writers) {
			try {
				writer.writeHeader();
			} catch (IOException e) {
				logger.error("Could not write header to file.", e);
			}
		}

		// Configure atomcache
		UserConfiguration cacheConfig = new UserConfiguration();
		if (pdbFilePath != null && !pdbFilePath.isEmpty()) {
			cacheConfig.setPdbFilePath(pdbFilePath);
			cacheConfig.setCacheFilePath(pdbFilePath);
		}
		AtomCache cache = new AtomCache(cacheConfig);
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
		cache.setFiletype(StructureFiletype.CIF);

		long startTime = System.nanoTime();

		// Start the workers in a fixed threaded pool
		ExecutorService executor = Executors.newFixedThreadPool(threads);
		for (String name : names) {
			StructureIdentifier id = new StructureName(name);
			Runnable worker = new QuatSymmWorker(id, sparams, cparams, cache,
					writers, show3d);
			executor.execute(worker);
		}
		executor.shutdown();
		while (!executor.isTerminated())
			Thread.sleep(100); // sleep .1 seconds

		long elapsed = (System.nanoTime() - startTime) / 1000000;
		long meanRT = (long) (elapsed / (float) names.size());
		logger.info("Total runtime: " + elapsed + ", mean runtime: " + meanRT);

		// Close any writers of output
		for (QuatSymmWriter writer : writers)
			writer.close();
	}

	/**
	 * Creates the options
	 *
	 * @param optionOrder
	 *            An empty map, which will be filled in with the option order
	 * @return all Options
	 */
	private static Options getOptions() {

		OptionGroup grp; // For mutually exclusive options
		Option opt;
		// Note: When adding an option, also add its long name to the

		Options options = new Options();
		options.addOption("h", "help", false, "Print usage information");
		options.addOption(Option.builder().longOpt("version").hasArg(false)
				.desc("Print CE-Symm version").build());

		// Input file
		options.addOption(Option.builder("i").longOpt("input").hasArg(true)
				.argName("file")
				.desc("File listing whitespace-delimited query structures")
				.build());

		// Logger control
		grp = new OptionGroup();
		opt = Option.builder("v").longOpt("verbose")
				.desc("Output verbose logging information.").build();
		grp.addOption(opt);
		opt = Option
				.builder("q")
				.longOpt("noverbose")
				.desc("Disable verbose logging information, as well as the "
						+ "default (--simple) output.")
				.build();
		grp.addOption(opt);
		options.addOptionGroup(grp);

		// Output formats
		options.addOption(Option
				.builder("o")
				.longOpt("stats")
				.hasArg()
				.optionalArg(true)
				.argName("file")
				.desc("Output a tsv file with detailed symmetry information (default)")
				.build());
		options.addOption(Option
				.builder("f")
				.longOpt("fasta")
				.hasArg()
				.optionalArg(true)
				.argName("file")
				.desc("Output alignment as FASTA alignment output")
				.build());

		// jmol
		grp = new OptionGroup();
		opt = Option
				.builder("j")
				.longOpt("show3d")
				.hasArg(false)
				.desc("Force jMol display for each structure "
						+ "[default for <10 structures when specified on command"
						+ " line]").build();
		grp.addOption(opt);
		opt = Option
				.builder("J")
				.longOpt("noshow3d")
				.hasArg(false)
				.desc("Disable jMol display [default with --input "
						+ "or for >=10 structures]").build();
		grp.addOption(opt);
		options.addOptionGroup(grp);

		// PDB_DIR
		options.addOption(Option
				.builder()
				.longOpt("pdbfilepath")
				.hasArg(true)
				.argName("dir")
				.desc("Download directory for new structures [default tmp folder]. "
						+ "Can also be set with the PDB_DIR environmental variable.")
				.build());

		options.addOption(Option.builder().longOpt("threads").hasArg(true)
				.desc("Number of threads [default cores-1]").build());

		// Parameters
		options.addOption(Option
				.builder()
				.longOpt("minSeqLen")
				.hasArg(true)
				.argName("int")
				.desc("The minimum subunit length to be considered for "
						+ "clustering and symmetry " + "analysis (default: 20)")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("minSeqId")
				.hasArg(true)
				.argName("float")
				.desc("Sequence identity threshold to consider for the sequence"
						+ " subunit clustering. Two subunits with sequence "
						+ "identity equal or higher than the threshold will be "
						+ "clustered together (range: [0,1], default: 0.95)")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("minSequenceCoverage")
				.hasArg(true)
				.argName("float")
				.desc("The minimum coverage of the sequence alignment between "
						+ "two subunits to be clustered together "
						+ "(range: [0,1], default: 0.9)").build());

		options.addOption(Option
				.builder()
				.longOpt("minStructureCoverage")
				.hasArg(true)
				.argName("float")
				.desc("The minimum coverage of the structure alignment between "
						+ "two subunits to be clustered together "
						+ "(range: [0,1], default: 0.9)").build());

		options.addOption(Option
				.builder()
				.longOpt("maxClustRmsd")
				.hasArg(true)
				.argName("float")
				.desc("Structure similarity threshold (measured with RMSD) "
						+ "to consider for the structural subunit clustering. "
						+ "(default: 3.0 A)").build());

		options.addOption(Option
				.builder()
				.longOpt("clustMethod")
				.hasArg(true)
				.argName("str")
				.desc("Method to cluster the subunits: "
						+ CliTools
								.getEnumValuesAsString(SubunitClustererMethod.class)
						+ "(default: STRUCTURE)").build());

		options.addOption(Option
				.builder()
				.longOpt("maxSymmRmsd")
				.hasArg(true)
				.argName("float")
				.desc("Structure similarity threshold (measured with RMSD) "
						+ "in the symmetry detection (default: 7.0 A)").build());

		return options;
	}

	/**
	 * Parse a whitespace-delimited file containing structure names
	 *
	 * @throws FileNotFoundException
	 */
	public static List<String> parseInputStructures(String filename)
			throws FileNotFoundException {
		File file = new File(filename);
		Scanner s = new Scanner(file);

		List<String> structures = new ArrayList<String>();
		while (s.hasNext()) {
			String name = s.next();
			if (name.startsWith("#")) {
				// comment
				s.nextLine();
			} else {
				structures.add(name);
			}
		}
		s.close();
		return structures;
	}

}
