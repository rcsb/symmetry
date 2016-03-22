package demo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.vecmath.Vector3d;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.ce.CeParameters.ScoringStrategy;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.CliTools;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.io.util.FileDownloadUtils;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.symmetry.gui.SymmetryDisplay;
import org.biojava.nbio.structure.symmetry.gui.SymmetryGui;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.OrderDetectorMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Main executable for running CE-Symm. Run with -h for usage help, or without
 * arguments for interactive mode using the {@link SymmetryGui}.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class CeSymmMain {

	private static final Logger logger = LoggerFactory
			.getLogger(CeSymmMain.class);

	public static void main(String[] args) throws InterruptedException {
		// Begin argument parsing
		final String usage = "[OPTIONS] [structures...]";
		final String header = "Determine the order for each structure, which may "
				+ "be PDB IDs, SCOP domains, or file paths. If none are given, the "
				+ "user will be prompted at startup.";
		final Map<String, Integer> optionOrder = new HashMap<String, Integer>();
		Options options = getOptions(optionOrder);
		CommandLineParser parser = new GnuParser();
		HelpFormatter help = new HelpFormatter();
		help.setOptionComparator(new Comparator<Option>() {
			@Override
			public int compare(Option o1, Option o2) {
				Integer i1 = optionOrder.get(o1.getLongOpt());
				Integer i2 = optionOrder.get(o2.getLongOpt());
				// Check that we didn't miss any options setting up optionOrder
				assert i1 != null : o1.getLongOpt();
				assert i2 != null : o2.getLongOpt();
				return i1.compareTo(i2);
			}
		});

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
			String version = CeSymmMain.class.getPackage()
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

		// Show jmol?
		// Default to false with --input or with >=10 structures
		boolean displayAlignment = !cli.hasOption("input") && names.size() < 10;
		if (cli.hasOption("noshow3d")) {
			displayAlignment = false;
		}
		if (cli.hasOption("show3d")) {
			displayAlignment = true;
		}

		// AtomCache options
		String pdbFilePath = null;
		if (cli.hasOption("pdbfilepath")) {
			pdbFilePath = cli.getOptionValue("pdbfilepath");
			pdbFilePath = FileDownloadUtils.expandUserHome(pdbFilePath);
		}

		// SCOP version
		if (cli.hasOption("scopversion")) {
			String scopVersion = cli.getOptionValue("scopversion");
			ScopFactory.setScopDatabase(scopVersion);
		}

		// Output formats
		List<CeSymmWriter> writers = new ArrayList<CeSymmWriter>();

		if (cli.hasOption("xml")) {
			String filename = cli.getOptionValue("xml");
			try {
				writers.add(new XMLWriter(filename));
			} catch (IOException e) {
				logger.error("Error: Ignoring file " + filename + ".");
				logger.error(e.getMessage());
			}
		}

		// Print summary to std out? verbose option
		// Default to false with --input or with >=10 structures
		boolean verbose = !cli.hasOption("input") && names.size() < 10;
		if (cli.hasOption("noverbose")) {
			verbose = false;
		}
		if (cli.hasOption("verbose")) {
			String filename = cli.getOptionValue("verbose");
			try {
				writers.add(new VerboseWriter(filename));
			} catch (IOException e) {
				logger.error("Error: Ignoring file " + filename + ".");
				logger.error(e.getMessage());
			}
		} else if (verbose) {
			try {
				writers.add(new VerboseWriter("-"));
			} catch (IOException e) {
				logger.error(e.getMessage());
			}
		}

		if (cli.hasOption("stats")) {
			String filename = cli.getOptionValue("stats");
			try {
				writers.add(new StatsWriter(filename));
			} catch (IOException e) {
				logger.error("Error: Ignoring file " + filename + ".");
				logger.error(e.getMessage());
			}
		}
		if (cli.hasOption("fatcat")) {
			String filename = cli.getOptionValue("fatcat");
			try {
				writers.add(new FatcatWriter(filename));
			} catch (IOException e) {
				logger.error("Error: Ignoring file " + filename + ".");
				logger.error(e.getMessage());
			}
		}
		if (cli.hasOption("fasta")) {
			String filename = cli.getOptionValue("fasta");
			try {
				writers.add(new FastaWriter(filename));
			} catch (IOException e) {
				logger.error("Error: Ignoring file " + filename + ".");
				logger.error(e.getMessage());
			}
		}

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
				logger.error("Illegal ordermethod. Requires on of "
						+ CliTools
								.getEnumValuesAsString(OrderDetectorMethod.class));
				System.exit(1);
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
		if (cli.hasOption("scorethreshold")) {
			String strVal = cli.getOptionValue("scorethreshold");
			try {
				double val = Double.parseDouble(strVal);
				if (val < 0) {
					logger.error("Invalid scorethreshold: " + strVal);
					System.exit(1);
				}
				params.setScoreThreshold(val);
			} catch (NumberFormatException e) {
				logger.error("Invalid scorethreshold: " + strVal);
				System.exit(1);
			}
		}
		if (cli.hasOption("ssethreshold")) {
			String strVal = cli.getOptionValue("ssethreshold");
			try {
				double val = Double.parseDouble(strVal);
				if (val < 0) {
					logger.error("Invalid ssethreshold: " + strVal);
					System.exit(1);
				}
				params.setScoreThreshold(val);
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

		// Done parsing arguments

		// Configure atomcache
		UserConfiguration cacheConfig = new UserConfiguration();
		if (pdbFilePath != null && !pdbFilePath.isEmpty()) {
			cacheConfig.setPdbFilePath(pdbFilePath);
			cacheConfig.setCacheFilePath(pdbFilePath);
		}
		AtomCache cache = new AtomCache(cacheConfig);
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);

		// Write the headers of the files
		for (CeSymmWriter writer : writers) {
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
			Runnable worker = new CeSymmWorker(id, params, cache, writers,
					displayAlignment);
			executor.execute(worker);
		}
		executor.shutdown();
		while (!executor.isTerminated())
			Thread.sleep(60000); // sleep 60 seconds

		long elapsed = (System.nanoTime() - startTime) / 1000000;
		long meanRT = (long) (elapsed / (float) names.size());
		logger.info("Total runtime: " + elapsed + ", mean runtime: " + meanRT);

		// Close any writers of output
		for (CeSymmWriter writer : writers)
			writer.close();
	}

	/**
	 * Creates the options
	 * 
	 * @param optionOrder
	 *            An empty map, which will be filled in with the option order
	 * @return all Options
	 */
	private static Options getOptions(Map<String, Integer> optionOrder) {

		OptionGroup grp;
		Option opt;
		// Note: When adding an option, also add its long name to the
		// optionOrder map
		int optionNum = 0;

		Options options = new Options();
		options.addOption("h", "help", false, "Print usage information");
		optionOrder.put("help", optionNum++);
		options.addOption(Option.builder().longOpt("version").hasArg(false)
				.desc("Print CE-Symm version").build());
		optionOrder.put("version", optionNum++);

		// Input file
		options.addOption(Option.builder("i")
				.longOpt("input")
				.hasArg(true)
				.argName("file")
				.desc(
						"File listing whitespace-delimited query structures")
				.build());
		optionOrder.put("input", optionNum++);
		// Output formats
		options.addOption(Option.builder()
				.longOpt("xml")
				.hasArg(true)
				.argName("file")
				.desc(
						"Output alignment as XML (use --xml=- for standard out).")
				.build());
		optionOrder.put("xml", optionNum++);
		options.addOption(Option.builder("v").longOpt("verbose").hasArg(true)
				.argName("file")
				.desc("Output verbose summary of the job results.")
				.build());
		optionOrder.put("verbose", optionNum++);
		options.addOption(Option.builder("o")
				.longOpt("stats")
				.hasArg(true)
				.argName("file")
				.desc(
						"Output a tsv file with the main symmetry info.")
				.build());
		// verbose ouput
		opt = Option.builder("q")
				.longOpt("noverbose")
				.hasArg(false)
				.desc(
						"Disable verbose summary to the std output."
								+ "[default for >=10 structures when specified on command"
								+ " line]").build();
		optionOrder.put(opt.getLongOpt(), optionNum++);
		//grp.addOption(opt);
		options.addOption(opt);
		optionOrder.put("stats", optionNum++);
		options.addOption(Option.builder().longOpt("fasta").hasArg(true)
				.argName("file")
				.desc("Output alignment as FASTA alignment output")
				.build());
		optionOrder.put("fasta", optionNum++);
		options.addOption(Option.builder().longOpt("fatcat").hasArg(true)
				.argName("file")
				.desc("Output alignment as FATCAT output").build());
		optionOrder.put("fatcat", optionNum++);

		// jmol
		grp = new OptionGroup();
		opt = Option.builder("j")
				.longOpt("show3d")
				.hasArg(false)
				.desc( "Force jMol display for each structure "
						+ "[default for <10 structures when specified on command"
						+ " line]")
				.build();
		grp.addOption(opt);
		optionOrder.put(opt.getLongOpt(), optionNum++);
		opt = Option.builder("J")
				.longOpt("noshow3d")
				.hasArg(false)
				.desc( "Disable jMol display [default with --input "
						+ "or for >=10 structures]")
				.build();
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		options.addOptionGroup(grp);

		// enums
		options.addOption(Option.builder()
				.longOpt("ordermethod")
				.hasArg(true)
				.argName("class")
				.desc(
						"Order detection method. Can be a "
								+ "full class name or a short class name from the "
								+ "org.biojava.nbio.structure.align.symmetry.internal package. "
								+ "[default SequenceFunctionOrderDetector]")
				.build());
		optionOrder.put("ordermethod", optionNum++);

		options.addOption(Option.builder()
				.longOpt("refinemethod")
				.hasArg(true)
				.argName("class")
				.desc(
						"Refiner method. Can be a "
								+ "full class name or a short class name from the "
								+ "org.biojava.nbio.structure.align.symmetry.internal package. "
								+ "[default Single]").build());
		optionOrder.put("refinemethod", optionNum++);

		options.addOption(Option.builder()
				.longOpt("symmtype")
				.hasArg(true)
				.argName("class")
				.desc(
						"Symmetry Type. Can be a "
								+ "full class name or a short class name from the "
								+ "org.biojava.nbio.structure.align.symmetry.internal package. "
								+ "[default Auto]").build());
		optionOrder.put("symmtype", optionNum++);

		// PDB_DIR
		options.addOption(Option.builder()
				.longOpt("pdbfilepath")
				.hasArg(true)
				.argName("dir")
				.desc(
						"Download directory for new "
								+ "structures. Equivalent to passing -DPDB_DIR=dir to the VM. "
								+ "[default temp folder]").build());
		optionOrder.put("pdbfilepath", optionNum++);


		grp = new OptionGroup();
		opt = Option.builder()
				.longOpt("pdbdirsplit")
				.hasArg(false)
				.desc(
						"Ignored. For backwards compatibility only. [default]")
				.build();
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		// grp.setSelected(opt);
		opt = Option.builder().longOpt("nopdbdirsplit").hasArg(false)
				.desc("Ignored. For backwards compatibility only.")
				.build();
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		options.addOptionGroup(grp);

		options.addOption(Option.builder().longOpt("threads").hasArg(true)
				.desc("Number of threads [default cores-1]")
				.build());
		optionOrder.put("threads", optionNum++);

		// Parameters
		options.addOption(Option.builder()
				.longOpt("maxgapsize")
				.hasArg(true)
				.argName("int")
				.desc(
						"This parameter configures the maximum gap size "
								+ "G, that is applied during the AFP extension. The "
								+ "larger the value, the longer the calculation time "
								+ "can become, Default value is 30. Set to 0 for no limit.")
				.build());
		optionOrder.put("maxgapsize", optionNum++);

		options.addOption(Option.builder()
				.longOpt("scoringstrategy")
				.hasArg(true)
				.argName("str")
				.desc(
						"Which scoring function to use: "
								+ CliTools
										.getEnumValuesAsString(ScoringStrategy.class))
				.build());
		optionOrder.put("scoringstrategy", optionNum++);

		options.addOption(Option.builder()
				.longOpt("winsize")
				.hasArg(true)
				.argName("int")
				.desc(
						"This configures the fragment size m of Aligned Fragment Pairs (AFPs).")
				.build());
		optionOrder.put("winsize", optionNum++);

		options.addOption(Option.builder()
				.longOpt("maxrmsd")
				.hasArg(true)
				.argName("float")
				.desc(
						"The maximum RMSD at which to stop alignment "
								+ "optimization. (default: unlimited=99)")
				.build());
		optionOrder.put("maxrmsd", optionNum++);

		options.addOption(Option.builder()
				.longOpt("gapopen")
				.hasArg(true)
				.argName("float")
				.desc(
						"Gap opening penalty during alignment optimization [default: 5.0].")
				.build());
		optionOrder.put("gapopen", optionNum++);

		options.addOption(Option.builder()
				.longOpt("gapextension")
				.hasArg(true)
				.argName("float")
				.desc(
						"Gap extension penalty during alignment optimization [default: 0.5].")
				.build());
		optionOrder.put("gapextension", optionNum++);

		options.addOption(Option.builder()
				.longOpt("symmlevels")
				.hasArg(true)
				.argName("int")
				.desc(
						"Run iteratively the algorithm to find multiple symmetry levels. "
								+ "This specifies the maximum number of symmetry levels. 0 means unbounded"
								+ " [default: 0].").build());
		optionOrder.put("symmlevels", optionNum++);

		options.addOption(Option.builder()
				.longOpt("noopt")
				.hasArg(false)
				.desc(
						"Disable optimization of the resulting symmetry alignment.")
				.build());
		optionOrder.put("noopt", optionNum++);

		options.addOption(Option.builder()
				.longOpt("scorethreshold")
				.hasArg(true)
				.argName("float")
				.desc(
						"The score threshold. TM-scores above this value "
								+ "will be considered significant results "
								+ "[default: 0.4, interval [0.0,1.0]].")
				.build());
		optionOrder.put("scorethreshold", optionNum++);

		options.addOption(Option.builder()
				.longOpt("ssethreshold")
				.hasArg(true)
				.argName("int")
				.desc(
						"The SSE threshold. Number of secondary structure"
								+ "elements for repeat below this value will be considered "
								+ "asymmetric results. 0 means unbounded."
								+ "[default: 0].").build());
		optionOrder.put("ssethreshold", optionNum++);

		options.addOption(Option.builder()
				.longOpt("maxorder")
				.hasArg(true)
				.argName("int")
				.desc(
						"The maximum number of symmetric repeats [default: 8].")
				.build());
		optionOrder.put("maxorder", optionNum++);

		options.addOption(Option.builder()
				.longOpt("rndseed")
				.hasArg(true)
				.argName("int")
				.desc(
						"The random seed used in optimization, for reproducibility "
								+ "of the results [default: 0].").build());
		optionOrder.put("rndseed", optionNum++);

		options.addOption(Option.builder()
				.longOpt("minlen")
				.hasArg(true)
				.argName("int")
				.desc(
						"The minimum length, expressed in number of core "
								+ "aligned residues, of a symmetric repeat [default: 15].")
				.build());
		optionOrder.put("minlen", optionNum++);

		options.addOption(Option.builder()
				.longOpt("dcutoff")
				.hasArg(true)
				.argName("float")
				.desc(
						"The maximum distance, in A, allowed between any two aligned "
								+ "residue positions [default: 7.0].")
				.build());
		optionOrder.put("dcutoff", optionNum++);

		options.addOption(Option.builder()
				.longOpt("scopversion")
				.hasArg(true)
				.argName("str")
				.desc(
						"Version of SCOP or SCOPe to use "
								+ "when resolving SCOP identifiers [defaults to latest SCOPe]")
				.build());
		optionOrder.put("scopversion", optionNum++);

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

	// Output formats

	/**
	 * Parent class for all output formats All methods are empty stubs, which
	 * should be overridden to write data to the writer.
	 */
	private static abstract class CeSymmWriter {
		protected PrintWriter writer;

		public CeSymmWriter(PrintWriter writer) {
			this.writer = writer;
		}

		public CeSymmWriter(String filename) throws IOException {
			this(openOutputFile(filename));
		}

		abstract public void writeHeader() throws IOException;

		abstract public void writeResult(CeSymmResult result)
				throws IOException;

		public void close() {
			if (writer != null) {
				writer.flush();
				writer.close();
			}
		}

		/**
		 * Opens 'filename' for writing.
		 * 
		 * @param filename
		 *            Name of output file, or '-' for standard out
		 * @throws IOException
		 */
		public static PrintWriter openOutputFile(String filename)
				throws IOException {
			if (filename.equals("-")) {
				return new PrintWriter(System.out, true);
			}
			return new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		}
	}

	private static class XMLWriter extends CeSymmWriter {
		public XMLWriter(String filename) throws IOException {
			super(filename);
		}

		@Override
		public void writeResult(CeSymmResult result) throws IOException {
			if (result.getMultipleAlignment() != null) {
				writer.append(MultipleAlignmentWriter.toXML(result
						.getMultipleAlignment().getEnsemble()));
				writer.flush();
			}
		}

		@Override
		public void writeHeader() throws IOException {
			// No header for XML file
		}
	}

	private static class FatcatWriter extends CeSymmWriter {
		public FatcatWriter(String filename) throws IOException {
			super(filename);
		}

		@Override
		public void writeResult(CeSymmResult result) {
			writer.write(MultipleAlignmentWriter.toFatCat(result
					.getMultipleAlignment()));
			writer.println("//");
			writer.flush();
		}

		@Override
		public void writeHeader() throws IOException {
			// No header for FatCat file
		}
	}

	private static class FastaWriter extends CeSymmWriter {
		public FastaWriter(String filename) throws IOException {
			super(filename);
		}

		@Override
		public void writeResult(CeSymmResult result) {
			writer.write(MultipleAlignmentWriter.toFASTA(result
					.getMultipleAlignment()));
			writer.println("//");
			writer.flush();
		}

		@Override
		public void writeHeader() throws IOException {
			// No header for Fasta files
		}
	}

	private static class StatsWriter extends CeSymmWriter {

		public StatsWriter(String filename) throws IOException {
			super(filename);
		}

		@Override
		public void writeHeader() {
			writer.println("Name\t" + "Order\t" + "SymmGroup\t" + "Refined\t"
					+ "Symm-Levels\t" + "Symm-Type\t" + "Rotation-Angle\t"
					+ "Screw-Translation\t" + "TMscore\t" + "RMSD\t"
					+ "Symm-TMscore\t" + "Symm-RMSD\t" + "RepeatLength\t"
					+ "Length\t" + "CoreLength\t" + "Coverage");
			writer.flush();
		}

		@Override
		public void writeResult(CeSymmResult result) throws IOException {

			String id = null;

			try {
				id = result.getStructureId().getIdentifier();
				int repeatLen = 0;
				int totalLen = result.getSelfAlignment().getOptLength();
				int coreLen = totalLen;
				double coverage = result.getSelfAlignment().getCoverage1() / 100;
				String group = result.getSymmGroup();
				int order = result.getSymmOrder();
				double symmrmsd = 0.0;
				double symmscore = 0.0;

				RotationAxis rot = new RotationAxis(result.getSelfAlignment()
						.getBlockRotationMatrix()[0], result.getSelfAlignment()
						.getBlockShiftVector()[0]);
				double rotation_angle = rot.getAngle() * 57.2957795;
				double screw_translation = new Vector3d(rot
						.getScrewTranslation().getCoords()).length();

				// If there is refinement alignment
				if (result.isRefined()) {
					MultipleAlignment msa = result.getMultipleAlignment();
					symmrmsd = msa.getScore(MultipleAlignmentScorer.RMSD);
					symmscore = msa
							.getScore(MultipleAlignmentScorer.AVGTM_SCORE);

					repeatLen = msa.length();
					double structureLen = result.getAtoms().length;
					totalLen = repeatLen * msa.size();
					coreLen = msa.getCoreLength() * msa.size();
					coverage = totalLen / structureLen;

					rot = new RotationAxis(result.getAxes().getElementaryAxes()
							.get(0));
					rotation_angle = rot.getAngle() * 57.2957795;
					screw_translation = new Vector3d(rot.getScrewTranslation()
							.getCoords()).length();
				}

				writer.format(
						"%s\t%d\t%s\t%b\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t"
								+ "%.2f\t%d\t%d\t%d\t%.2f\n", id, order, group,
						result.isRefined(), result.getSymmLevels(), result
								.getType(), rotation_angle, screw_translation,
						result.getSelfAlignment().getTMScore(), result
								.getSelfAlignment().getTotalRmsdOpt(),
						symmscore, symmrmsd, repeatLen, totalLen, coreLen,
						coverage);
			} catch (Exception e) {
				// If any exception occurs when writing the results store empty
				// better
				logger.warn("Could not write result... storing empty row.", e);
				writeEmptyRow(writer, id);
			}

			writer.flush();
		}

		private void writeEmptyRow(PrintWriter writer, String id) {
			writer.format(
					"%s\t%d\t%s\t%b\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t"
							+ "%.2f\t%d\t%d\t%d\t%.2f\n", id, 1, "C1", false,
					0, SymmetryType.DEFAULT, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
					0, 0, 0.0);
		}
	}

	private static class VerboseWriter extends CeSymmWriter {

		public VerboseWriter(String filename) throws IOException {
			super(filename);
		}

		@Override
		public void writeHeader() {
			writer.println("Summary of the job results:");
			writer.flush();
		}

		@Override
		public void writeResult(CeSymmResult result) throws IOException {

			writer.append(result.getStructureId().getIdentifier());

			if (result.isRefined()) {
				writer.append(" could be refined into symmetry order "
						+ result.getMultipleAlignment().size()
						+ ". The structure is ");
				if (result.isSignificant()) {
					writer.append("significant (symmetric).");
				} else
					writer.append("not significant (asymmetric).");
				writer.append(System.getProperty("line.separator"));
			} else {
				writer.append(" could not be refined. "
						+ "The structure has no symmetry (asymmetric)");
				writer.append(System.getProperty("line.separator"));
			}
			writer.flush();
		}
	}

	/**
	 * This Runnable implementation runs CeSymm on the input structure and with
	 * the input parameters and writes the results to the output writers. If the
	 * 3D visualization is turned on, it creates a new thread with the Jmol
	 * frame.
	 * 
	 * @author Aleix Lafita
	 *
	 */
	private static class CeSymmWorker implements Runnable {

		private StructureIdentifier id;
		private CESymmParameters params;
		private AtomCache cache;
		private List<CeSymmWriter> writers;
		private boolean show3d;

		public CeSymmWorker(StructureIdentifier id, CESymmParameters params,
				AtomCache cache, List<CeSymmWriter> writers, boolean show3d) {
			this.id = id;
			this.cache = cache;
			this.writers = writers;
			this.show3d = show3d;
			this.params = params;
		}

		@Override
		public void run() {

			try {
				// Obtain the structure representation
				Structure s = cache.getStructure(id);
				Atom[] atoms = StructureTools.getRepresentativeAtomArray(s);

				// Run the symmetry analysis
				CeSymmResult result = CeSymm.analyze(atoms, params);

				// Write into the output files
				for (CeSymmWriter writer : writers) {
					try {
						synchronized (writer) {
							writer.writeResult(result);
						}
					} catch (Exception e) {
						logger.error(
								"Could not save results for "
										+ id.getIdentifier(), e);
					}
				}

				// Display alignment in 3D Jmol
				if (show3d) {
					try {
						SymmetryDisplay.display(result);
					} catch (Exception e) {
						//TODO NullPointer when changing title fixed in biojava 5.0
					}
				}
			} catch (IOException e) {
				logger.error("Could not load Structure " + id.getIdentifier(),
						e);
			} catch (Exception e) {
				logger.error("Could not complete job: " + id.getIdentifier(), e);
			}
			logger.info("Finished job: " + id);
		}
	}

}
