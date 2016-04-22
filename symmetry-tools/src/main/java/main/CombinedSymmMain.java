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
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.align.ce.CeParameters.ScoringStrategy;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.align.util.CliTools;
import org.biojava.nbio.structure.io.util.FileDownloadUtils;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.gui.SymmetryGui;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.OrderDetectorMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import workers.CombinedSymmWorker;
import writers.CombinedSymmWriter;

/**
 * This script calculates the prevalence of internal and quaternary symmetry
 * combination across a list of input entries.
 * <p>
 * Internal and quaternary symmetry combination can be found in: (1) uneven
 * stoichiometries where the lower stoichiometry chain is internally symmetric
 * and the complex is pseudosymmetric after internal symmetry split, and (2)
 * point group symmetry that is amplified by internal symmetry.
 * <p>
 * The script downloads the biological assembly 1 from the PDB, runs the
 * {@link QuatSymmetryDetector} algorithm first and then {@link CeSymm} on each
 * chain of the assembly. Chains are split according to their internal symmetry
 * and the {@link QuatSymmetryDetector} is run again on the complex.
 * <p>
 * The following outcomes are possible: (1) The protein was monomeric and
 * internally symmetric, (2) the complex was heteromeric and asymmetric and
 * gained overall symmetry with internal symmetry, uneven stoichiometry, and (3)
 * the complex was homomeric and symmetric and internal symmetry amplified the
 * point group order.
 * 
 * @author Aleix Lafita
 * @since Apr 2016
 * @version beta testing
 *
 */
public class CombinedSymmMain {

	private static final Logger logger = LoggerFactory
			.getLogger(CombinedSymmMain.class);

	public static void main(String[] args) throws InterruptedException {
		// Begin argument parsing
		final String usage = "[OPTIONS] [structures...]";
		final String header = "Determine the order for each structure, which may "
				+ "be PDB IDs, SCOP domains, or file paths. If none are given, the "
				+ "user will be prompted at startup.";
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
			String version = CombinedSymmMain.class.getPackage()
					.getImplementationVersion();
			if (version == null || version.isEmpty()) {
				version = "(custom version)";
			}
			System.out.println("CE-Symm " + version);
			System.exit(0);
			return;
		}

		// input structures
		List<String> names;
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
				// No structures given; prompt user with GUI
				SymmetryGui.getInstance();
				return;
			} else {
				// take names from the command line arguments
				names = Arrays.asList(args);
			}
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
		List<CombinedSymmWriter> writers = new ArrayList<CombinedSymmWriter>();

		if (cli.hasOption("output")) {
			String filename = cli.getOptionValue("output");
			if (filename == null || filename.isEmpty())
				filename = "-"; // standard out
			try {
				writers.add(new CombinedSymmWriter(filename));
			} catch (IOException e) {
				logger.error("Error: Ignoring file " + filename + ".");
				logger.error(e.getMessage());
			}
		}

		// Default Writer
		if (writers.isEmpty() && !cli.hasOption("noverbose")) {
			try {
				writers.add(new CombinedSymmWriter("-"));
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

		CESymmParameters params = new CESymmParameters();

		if (cli.hasOption("maxgapsize")) {
			String gapStr = cli.getOptionValue("maxgapsize");
			try {
				int gap = Integer.parseInt(gapStr);
				if (gap < 1) {
					logger.error("Invalid maxgapsize: " + gap);
					System.exit(1);
				}
				params.setMaxGapSize(gap);
			} catch (NumberFormatException e) {
				logger.error("Invalid maxgapsize: " + gapStr);
				System.exit(1);
			}
		}
		if (cli.hasOption("scoringstrategy")) {
			String stratStr = cli.getOptionValue("scoringstrategy");
			ScoringStrategy strat;
			try {
				strat = ScoringStrategy.valueOf(stratStr.toUpperCase());
				params.setScoringStrategy(strat);
			} catch (IllegalArgumentException e) {
				// give up
				logger.error("Illegal scoringstrategy. Requires on of "
						+ CliTools.getEnumValuesAsString(ScoringStrategy.class));
				System.exit(1);
			}
		}
		if (cli.hasOption("winsize")) {
			String winStr = cli.getOptionValue("winsize");
			try {
				int win = Integer.parseInt(winStr);
				if (win < 1) {
					logger.error("Invalid winsize: " + winStr);
					System.exit(1);
				}
				params.setWinSize(win);
			} catch (NumberFormatException e) {
				logger.error("Invalid winsize: " + winStr);
				System.exit(1);
			}

		}
		if (cli.hasOption("maxrmsd")) {
			String strVal = cli.getOptionValue("maxrmsd");
			try {
				double val = Double.parseDouble(strVal);
				if (val < 0) {
					logger.error("Invalid maxrmsd: " + strVal);
					System.exit(1);
				}
				params.setMaxOptRMSD(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid maxrmsd: " + strVal);
				System.exit(1);
			}

		}

		if (cli.hasOption("nointernalgaps")) {
			params.setGaps(false);
		}
		if (cli.hasOption("gapopen")) {
			String strVal = cli.getOptionValue("gapopen");
			try {
				double val = Double.parseDouble(strVal);
				if (val < 0) {
					logger.error("Invalid gapopen: " + strVal);
					System.exit(1);
				}
				params.setGapOpen(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid gapopen: " + strVal);
				System.exit(1);
			}

		}
		if (cli.hasOption("gapextension")) {
			String strVal = cli.getOptionValue("gapextension");
			try {
				double val = Double.parseDouble(strVal);
				if (val < 0) {
					logger.error("Invalid gapextension: " + strVal);
					System.exit(1);
				}
				params.setGapExtension(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid gapextension: " + strVal);
				System.exit(1);
			}
		}
		if (cli.hasOption("ordermethod")) {
			String strVal = cli.getOptionValue("ordermethod");
			OrderDetectorMethod val;
			try {
				val = OrderDetectorMethod.valueOf(strVal.toUpperCase());
				params.setOrderDetectorMethod(val);
			} catch (IllegalArgumentException e) {
				// give up
				logger.error("Illegal ordermethod. Requires on of {}", CliTools
						.getEnumValuesAsString(OrderDetectorMethod.class));
				System.exit(1);
			}
		}
		if (cli.hasOption("order")) {
			String strVal = cli.getOptionValue("order");
			try {
				int val = Integer.parseInt(strVal);
				params.setUserOrder(val);
				if (val > 0) {
					params.setOrderDetectorMethod(OrderDetectorMethod.USER_INPUT);
				}
			} catch (NumberFormatException e) {
				logger.error("Invalid order: " + strVal);
			}
		}
		if (cli.hasOption("refinemethod")) {
			String strVal = cli.getOptionValue("refinemethod");
			RefineMethod val;
			try {
				val = RefineMethod.valueOf(strVal.toUpperCase());
				params.setRefineMethod(val);
			} catch (IllegalArgumentException e) {
				// give up
				logger.error("Illegal refinemethod. Requires on of "
						+ CliTools.getEnumValuesAsString(RefineMethod.class));
				System.exit(1);
			}
		}
		if (cli.hasOption("symmtype")) {
			String strVal = cli.getOptionValue("symmtype");
			SymmetryType val;
			try {
				val = SymmetryType.valueOf(strVal.toUpperCase());
				params.setSymmType(val);
			} catch (IllegalArgumentException e) {
				// give up
				logger.error("Illegal symmtype. Requires on of "
						+ CliTools.getEnumValuesAsString(RefineMethod.class));
				System.exit(1);
			}
		}
		if (cli.hasOption("maxorder")) {
			String strVal = cli.getOptionValue("maxorder");
			try {
				int val = Integer.parseInt(strVal);
				if (val < 0) {
					logger.error("Invalid maxorder: " + strVal);
					System.exit(1);
				}
				params.setMaxSymmOrder(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid maxorder: " + strVal);
				System.exit(1);
			}
		}
		if (cli.hasOption("rndseed")) {
			String strVal = cli.getOptionValue("rndseed");
			try {
				int val = Integer.parseInt(strVal);
				if (val < 0) {
					logger.error("Invalid rndseed: " + strVal);
					System.exit(1);
				}
				params.setRndSeed(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid rndseed: " + strVal);
				System.exit(1);
			}
		}
		boolean optimize = cli.hasOption("opt") || !cli.hasOption("noopt");
		params.setOptimization(optimize);
		if (cli.hasOption("symmlevels")) {
			String strVal = cli.getOptionValue("symmlevels");
			try {
				int val = Integer.parseInt(strVal);
				params.setSymmLevels(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid symmlevels: " + strVal);
				System.exit(1);
			}
		}
		if (cli.hasOption("unrefinedscorethreshold")) {
			String strVal = cli.getOptionValue("unrefinedscorethreshold");
			try {
				double val = Double.parseDouble(strVal);
				if (val < 0) {
					logger.error("Invalid unrefinedscorethreshold: " + strVal);
					System.exit(1);
				}
				params.setUnrefinedScoreThreshold(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid unrefinedscorethreshold: " + strVal);
				System.exit(1);
			}
		}
		if (cli.hasOption("refinedscorethreshold")) {
			String strVal = cli.getOptionValue("refinedscorethreshold");
			try {
				double val = Double.parseDouble(strVal);
				if (val < 0) {
					logger.error("Invalid refinedscorethreshold: " + strVal);
					System.exit(1);
				}
				params.setRefinedScoreThreshold(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid refinedscorethreshold: " + strVal);
				System.exit(1);
			}
		}
		if (cli.hasOption("ssethreshold")) {
			String strVal = cli.getOptionValue("ssethreshold");
			try {
				int val = Integer.parseInt(strVal);
				if (val < 0) {
					logger.error("Invalid ssethreshold: " + strVal);
					System.exit(1);
				}
				params.setSSEThreshold(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid ssethreshold: " + strVal);
				System.exit(1);
			}
		}
		if (cli.hasOption("dcutoff")) {
			String strVal = cli.getOptionValue("dcutoff");
			try {
				double val = Double.parseDouble(strVal);
				if (val < 0) {
					logger.error("Invalid dcutoff: " + strVal);
					System.exit(1);
				}
				params.setDistanceCutoff(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid dcutoff: " + strVal);
				System.exit(1);
			}
		}
		if (cli.hasOption("minlen")) {
			String strVal = cli.getOptionValue("minlen");
			try {
				int val = Integer.parseInt(strVal);
				if (val < 0) {
					logger.error("Invalid minlen: " + strVal);
					System.exit(1);
				}
				params.setMinCoreLength(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid minlen: " + strVal);
				System.exit(1);
			}
		}

		verifyParams(params);

		QuatSymmetryParameters qparams = new QuatSymmetryParameters();

		// Done parsing arguments

		// Write the headers of the files
		for (CombinedSymmWriter writer : writers) {
			try {
				writer.writeHeader();
			} catch (IOException e) {
				logger.error("Could not write header to file.", e);
			}
		}
		long startTime = System.nanoTime();

		// Start the workers in a fixed threaded pool
		ExecutorService executor = Executors.newFixedThreadPool(threads);
		for (String name : names) {
			StructureIdentifier id = new StructureName(name);
			Runnable worker = new CombinedSymmWorker(id, params, qparams,
					writers);
			executor.execute(worker);
		}
		executor.shutdown();
		while (!executor.isTerminated())
			Thread.sleep(100); // sleep .1 seconds

		long elapsed = (System.nanoTime() - startTime) / 1000000;
		long meanRT = (long) (elapsed / (float) names.size());
		logger.info("Total runtime: " + elapsed + ", mean runtime: " + meanRT);

		// Close any writers of output
		for (CombinedSymmWriter writer : writers)
			writer.close();
	}

	/**
	 * Check dependencies between parameters
	 * 
	 * @param params
	 */
	private static void verifyParams(CESymmParameters params) {
		int order = params.getUserOrder();
		OrderDetectorMethod orderdetector = params.getOrderDetectorMethod();
		if (order > 0 && orderdetector != OrderDetectorMethod.USER_INPUT) {
			logger.info("--order={} is incompatible with --orderdetector={}",
					order, orderdetector);
			System.exit(1);
		}
		if (order < 1 && orderdetector == OrderDetectorMethod.USER_INPUT) {
			logger.info("USER_INPUT detector requires --order", orderdetector);
			System.exit(1);
		}
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
		// .hasArg(true)
		// .argName("file")
				.desc("Output verbose logging information.").build();
		grp.addOption(opt);
		opt = Option
				.builder("q")
				.longOpt("noverbose")
				.hasArg(false)
				.desc("Disable verbose logging information, as well as the default (--simple) output.")
				.build();
		grp.addOption(opt);
		options.addOptionGroup(grp);

		// Output formats
		options.addOption(Option.builder("o").longOpt("simple").hasArg()
				.optionalArg(true).argName("file")
				.desc("Output result in a simple format (default).").build());
		options.addOption(Option.builder().longOpt("stats").hasArg()
				.optionalArg(true).argName("file")
				.desc("Output a tsv file with detailed symmetry info.").build());
		options.addOption(Option
				.builder()
				.longOpt("tsv")
				.hasArg()
				.optionalArg(true)
				.argName("file")
				.desc("Output alignment as a tsv-formated list of aligned residues.")
				.build());
		options.addOption(Option
				.builder()
				.longOpt("xml")
				.hasArg()
				.optionalArg(true)
				.argName("file")
				.desc("Output alignment as XML (use --xml=- for standard out).")
				.build());
		options.addOption(Option.builder().longOpt("fatcat").hasArg()
				.optionalArg(true).argName("file")
				.desc("Output alignment as FATCAT output").build());
		options.addOption(Option.builder().longOpt("fasta").hasArg()
				.optionalArg(true).argName("file")
				.desc("Output alignment as FASTA alignment output").build());

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

		// enums
		options.addOption(Option
				.builder()
				.longOpt("ordermethod")
				.hasArg(true)
				.argName("str")
				.desc("Order detection method: "
						+ CliTools
								.getEnumValuesAsString(OrderDetectorMethod.class))
				.build());
		options.addOption(Option
				.builder()
				.longOpt("order")
				.hasArg(true)
				.argName("int")
				.desc("Force a particular order. If positive, implies --ordermethod=USER_INPUT.")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("refinemethod")
				.hasArg(true)
				.argName("str")
				.desc("Refiner method: "
						+ CliTools.getEnumValuesAsString(RefineMethod.class))
				.build());

		options.addOption(Option
				.builder()
				.longOpt("symmtype")
				.hasArg(true)
				.argName("str")
				.desc("Symmetry Type: "
						+ CliTools.getEnumValuesAsString(SymmetryType.class))
				.build());

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
				.longOpt("maxgapsize")
				.hasArg(true)
				.argName("int")
				.desc("This parameter configures the maximum gap size "
						+ "G, that is applied during the AFP extension. The "
						+ "larger the value, the longer the calculation time "
						+ "can become, Default value is 30. Set to 0 for no limit.")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("scoringstrategy")
				.hasArg(true)
				.argName("str")
				.desc("Which scoring function to use: "
						+ CliTools.getEnumValuesAsString(ScoringStrategy.class))
				.build());

		options.addOption(Option
				.builder()
				.longOpt("winsize")
				.hasArg(true)
				.argName("int")
				.desc("This configures the fragment size m of Aligned Fragment Pairs (AFPs).")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("maxrmsd")
				.hasArg(true)
				.argName("float")
				.desc("The maximum RMSD at which to stop alignment "
						+ "optimization. (default: unlimited=99)").build());

		options.addOption(Option
				.builder()
				.longOpt("nointernalgaps")
				.desc("Force alignment to include a residue from all repeats. "
						+ "(By default only 50% of repeats must be aligned in each column.)")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("gapopen")
				.hasArg(true)
				.argName("float")
				.desc("Gap opening penalty during alignment optimization [default: 5.0].")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("gapextension")
				.hasArg(true)
				.argName("float")
				.desc("Gap extension penalty during alignment optimization [default: 0.5].")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("symmlevels")
				.hasArg(true)
				.argName("int")
				.desc("Run iteratively the algorithm to find multiple symmetry levels. "
						+ "This specifies the maximum number of symmetry levels. 0 means unbounded"
						+ " [default: 0].").build());

		options.addOption(Option
				.builder()
				.longOpt("noopt")
				.hasArg(false)
				.desc("Disable optimization of the resulting symmetry alignment.")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("unrefinedscorethreshold")
				.hasArg(true)
				.argName("float")
				.desc("The unrefined score threshold. TM-scores above this "
						+ "value in the optimal self-alignment of the "
						+ "structure (before refinement) "
						+ "will be considered significant results "
						+ "[default: 0.4, interval [0.0,1.0]].").build());

		options.addOption(Option
				.builder()
				.longOpt("refinedscorethreshold")
				.hasArg(true)
				.argName("float")
				.desc("The refined score threshold. TM-scores above this "
						+ "value in the refined multiple alignment "
						+ "of repeats (after refinement) "
						+ "will be considered significant results "
						+ "[default: 0.36, interval [0.0,1.0]].").build());

		options.addOption(Option
				.builder()
				.longOpt("ssethreshold")
				.hasArg(true)
				.argName("int")
				.desc("The SSE threshold. Number of secondary structure "
						+ "elements for repeat below this value will be considered "
						+ "asymmetric results. 0 means unbounded."
						+ "[default: 0].").build());

		options.addOption(Option.builder().longOpt("maxorder").hasArg(true)
				.argName("int")
				.desc("The maximum number of symmetric repeats [default: 8].")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("rndseed")
				.hasArg(true)
				.argName("int")
				.desc("The random seed used in optimization, for reproducibility "
						+ "of the results [default: 0].").build());

		options.addOption(Option
				.builder()
				.longOpt("minlen")
				.hasArg(true)
				.argName("int")
				.desc("The minimum length, expressed in number of core "
						+ "aligned residues, of a symmetric repeat [default: 15].")
				.build());

		options.addOption(Option
				.builder()
				.longOpt("dcutoff")
				.hasArg(true)
				.argName("float")
				.desc("The maximum distance, in A, allowed between any two aligned "
						+ "residue positions [default: 7.0].").build());

		options.addOption(Option
				.builder()
				.longOpt("scopversion")
				.hasArg(true)
				.argName("str")
				.desc("Version of SCOP or SCOPe to use "
						+ "when resolving SCOP identifiers [defaults to latest SCOPe]")
				.build());

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
