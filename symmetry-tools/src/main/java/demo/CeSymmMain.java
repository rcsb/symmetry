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
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.swing.JOptionPane;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.ce.CeParameters.ScoringStrategy;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.CliTools;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.util.FileDownloadUtils;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.symmetry.gui.SymmetryDisplay;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.OrderDetectorMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.SymmetryAxes;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Main executable for running CE-Symm
 * Run with -h for usage, or without arguments for interactive mode
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class CeSymmMain {

	private static final Logger logger = LoggerFactory.getLogger(CeSymmMain.class);

	public static void main(String[] args) {
		// Begin argument parsing
		final String usage = "[OPTIONS] [structures...]";
		final String header = "Determine the order for each structure, which may " +
				"be PDB IDs, SCOP domains, or file paths. If none are given, the " +
				"user will be prompted at startup.";
		final Map<String,Integer> optionOrder = new HashMap<String,Integer>();
		Options options = getOptions(optionOrder);
		CommandLineParser parser = new GnuParser();
		HelpFormatter help = new HelpFormatter();
		help.setOptionComparator(new Comparator<Option>() {
			@Override
			public int compare(Option o1, Option o2) {
				Integer i1 = optionOrder.get(o1.getLongOpt());
				Integer i2 = optionOrder.get(o2.getLongOpt());
				// Check that we didn't miss any options setting up optionOrder
				assert i1!=null : o1.getLongOpt();
				assert i2!=null : o2.getLongOpt();
				return i1.compareTo(i2);
			}
		} );
		
		final CommandLine cli;
		try {
			cli = parser.parse(options,args,false);
		} catch (ParseException e) {
			System.err.println("Error: "+e.getMessage());
			help.printHelp(usage, header, options, "");
			System.exit(1);
			return;
		}

		args = cli.getArgs();

		// help
		if(cli.hasOption("help")) {
			help.printHelp(usage, header, options, "");
			System.exit(0);
			return;
		}
		// version
		if(cli.hasOption("version")) {
			String version = CeSymmMain.class.getPackage().getImplementationVersion();
			if(version == null || version.isEmpty()) {
				version = "(custom version)";
			}
			System.out.println("CE-Symm "+version);
			System.exit(0);
			return;
		}

		// input structures
		List<String> names;
		if(cli.hasOption("input")) {
			// read from file
			try {
				names = parseInputStructures(cli.getOptionValue("input"));
			} catch (FileNotFoundException e) {
				System.err.println("Error: File not found: "+cli.getOptionValue("input"));
				System.exit(1);
				return;
			}
			// append cli arguments
			names.addAll(Arrays.asList(args));
		} else {
			if(args.length == 0) {
				// No structures given; prompt user
				String name = promptUserForStructure();
				if( name == null) {
					//cancel
					return;
				}
				names = Arrays.asList(name);
			} else {
				// take names from the command line arguments
				names = Arrays.asList(args);
			}
		}

		// Show jmol?
		// Default to false with --input or with >=10 structures
		boolean displayAlignment = !cli.hasOption("input") && names.size() < 10;
		if(cli.hasOption("noshow3d")) {
			displayAlignment = false;
		}
		if(cli.hasOption("show3d")) {
			displayAlignment = true;
		}

		// AtomCache options
		String pdbFilePath = null;
		if( cli.hasOption("pdbfilepath") ) {
			pdbFilePath = cli.getOptionValue("pdbfilepath");
			pdbFilePath = FileDownloadUtils.expandUserHome(pdbFilePath);
		}

		// SCOP version
		if( cli.hasOption("scopversion")) {
			String scopVersion = cli.getOptionValue("scopversion");
			ScopFactory.setScopDatabase(scopVersion);
		}

		// Output formats
		List<CeSymmWriter> writers = new ArrayList<CeSymmWriter>();

		if(cli.hasOption("xml")) {
			String filename = cli.getOptionValue("xml");
			try {
				writers.add(new XMLWriter(filename));
			} catch (IOException e) {
				System.err.println("Error: Ignoring file "+filename+".");
				System.err.println(e.getMessage());
			}
		}
		if(cli.hasOption("stats")) {
			String filename = cli.getOptionValue("stats");
			try {
				writers.add(new StatsWriter(filename));
			} catch (IOException e) {
				System.err.println("Error: Ignoring file "+filename+".");
				System.err.println(e.getMessage());
			}
		}
		if(cli.hasOption("fatcat")) {
			String filename = cli.getOptionValue("fatcat");
			try {
				writers.add(new FatcatWriter(filename));
			} catch (IOException e) {
				System.err.println("Error: Ignoring file "+filename+".");
				System.err.println(e.getMessage());
			}
		}
		if(cli.hasOption("fasta")) {
			String filename = cli.getOptionValue("fasta");
			try {
				writers.add(new FastaWriter(filename));
			} catch (IOException e) {
				System.err.println("Error: Ignoring file "+filename+".");
				System.err.println(e.getMessage());
			}
		}

		Integer threads = 1;
		if(cli.hasOption("threads")) {
			threads = new Integer(cli.getOptionValue("threads"));
		}

		CESymmParameters params = new CESymmParameters();

		if(cli.hasOption("maxgapsize")) {
			String gapStr = cli.getOptionValue("maxgapsize");
			try {
				int gap = Integer.parseInt(gapStr);
				if(gap < 1) {
					System.err.println("Invalid maxgapsize: "+gap);
					System.exit(1);
				}
				params.setMaxGapSize(gap);
			} catch( NumberFormatException e) {
				System.err.println("Invalid maxgapsize: "+gapStr);
				System.exit(1);
			}
		}
		if(cli.hasOption("scoringstrategy")) {
			String stratStr = cli.getOptionValue("scoringstrategy");
			ScoringStrategy strat;
			try {
				strat = ScoringStrategy.valueOf( stratStr.toUpperCase());
				params.setScoringStrategy(strat);
			} catch (IllegalArgumentException e) {
				//give up
				System.err.println("Illegal scoringstrategy. Requires on of "+
						CliTools.getEnumValuesAsString(ScoringStrategy.class));
				System.exit(1);
			}
		}
		if(cli.hasOption("winsize")) {
			String winStr = cli.getOptionValue("winsize");
			try {
				int win = Integer.parseInt(winStr);
				if(win < 1) {
					System.err.println("Invalid winsize: "+winStr);
					System.exit(1);
				}
				params.setWinSize(win);
			} catch( NumberFormatException e) {
				System.err.println("Invalid winsize: "+winStr);
				System.exit(1);
			}

		}
		if(cli.hasOption("maxrmsd")) {
			String strVal = cli.getOptionValue("maxrmsd");
			try {
				double val = Double.parseDouble(strVal);
				if(val < 0) {
					System.err.println("Invalid maxrmsd: "+strVal);
					System.exit(1);
				}
				params.setMaxOptRMSD(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid maxrmsd: "+strVal);
				System.exit(1);
			}

		}
		if(cli.hasOption("gapopen")) {
			String strVal = cli.getOptionValue("gapopen");
			try {
				double val = Double.parseDouble(strVal);
				if(val < 0) {
					System.err.println("Invalid gapopen: "+strVal);
					System.exit(1);
				}
				params.setGapOpen(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid gapopen: "+strVal);
				System.exit(1);
			}

		}
		if(cli.hasOption("gapextension")) {
			String strVal = cli.getOptionValue("gapextension");
			try {
				double val = Double.parseDouble(strVal);
				if(val < 0) {
					System.err.println("Invalid gapextension: "+strVal);
					System.exit(1);
				}
				params.setGapExtension(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid gapextension: "+strVal);
				System.exit(1);
			}
		}
		if(cli.hasOption("ordermethod")) {
			String strVal = cli.getOptionValue("ordermethod");
			OrderDetectorMethod val;
			try {
				val = OrderDetectorMethod.valueOf( strVal.toUpperCase());
				params.setOrderDetectorMethod(val);
			} catch (IllegalArgumentException e) {
				//give up
				System.err.println("Illegal ordermethod. Requires on of "+
						CliTools.getEnumValuesAsString(OrderDetectorMethod.class));
				System.exit(1);
			}
		}
		if(cli.hasOption("refinemethod")) {
			String strVal = cli.getOptionValue("refineresult");
			RefineMethod val;
			try {
				val = RefineMethod.valueOf(strVal.toUpperCase());
				params.setRefineMethod(val);
			} catch (IllegalArgumentException e) {
				//give up
				System.err.println("Illegal refinemethod. Requires on of "+
						CliTools.getEnumValuesAsString(RefineMethod.class));
				System.exit(1);
			}
		}
		if(cli.hasOption("symmtype")) {
			String strVal = cli.getOptionValue("symmtype");
			SymmetryType val;
			try {
				val = SymmetryType.valueOf(strVal.toUpperCase());
				params.setSymmetryType(val);
			} catch (IllegalArgumentException e) {
				//give up
				System.err.println("Illegal symmtype. Requires on of "+
						CliTools.getEnumValuesAsString(RefineMethod.class));
				System.exit(1);
			}
		}
		if(cli.hasOption("maxorder")) {
			String strVal = cli.getOptionValue("maxorder");
			try {
				int val = Integer.parseInt(strVal);
				if(val < 0) {
					System.err.println("Invalid maxorder: "+strVal);
					System.exit(1);
				}
				params.setMaxSymmOrder(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid maxorder: "+strVal);
				System.exit(1);
			}
		}
		if(cli.hasOption("rndseed")) {
			String strVal = cli.getOptionValue("rndseed");
			try {
				int val = Integer.parseInt(strVal);
				if(val < 0) {
					System.err.println("Invalid rndseed: "+strVal);
					System.exit(1);
				}
				params.setSeed(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid rndseed: "+strVal);
				System.exit(1);
			}
		}
		if(cli.hasOption("opt")) {
			String strVal = cli.getOptionValue("opt");
			try {
				boolean val = Boolean.parseBoolean(strVal);
				params.setOptimization(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid opt: "+strVal);
				System.exit(1);
			}
		}
		if(cli.hasOption("multaxes")) {
			String strVal = cli.getOptionValue("multaxes");
			try {
				boolean val = Boolean.parseBoolean(strVal);
				params.setMultipleAxes(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid multaxes: "+strVal);
				System.exit(1);
			}
		}
		if(cli.hasOption("threshold")) {
			String strVal = cli.getOptionValue("threshold");
			try {
				double val = Double.parseDouble(strVal);
				if(val < 0) {
					System.err.println("Invalid threshold: "+strVal);
					System.exit(1);
				}
				params.setSymmetryThreshold(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid threshold: "+strVal);
				System.exit(1);
			}
		}
		if(cli.hasOption("dcutoff")) {
			String strVal = cli.getOptionValue("dcutoff");
			try {
				double val = Double.parseDouble(strVal);
				if(val < 0) {
					System.err.println("Invalid dcutoff: "+strVal);
					System.exit(1);
				}
				params.setDistanceCutoff(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid dcutoff: "+strVal);
				System.exit(1);
			}
		}
		if(cli.hasOption("minlen")) {
			String strVal = cli.getOptionValue("minlen");
			try {
				int val = Integer.parseInt(strVal);
				if(val < 0) {
					System.err.println("Invalid minlen: "+strVal);
					System.exit(1);
				}
				params.setMinSubunitLength(val);
			} catch( NumberFormatException e) {
				System.err.println("Invalid minlen: "+strVal);
				System.exit(1);
			}
		}

		// Done parsing arguments

		// Configure atomcache
		UserConfiguration cacheConfig = new UserConfiguration();
		if(pdbFilePath != null && !pdbFilePath.isEmpty()) {
			cacheConfig.setPdbFilePath(pdbFilePath);
			cacheConfig.setCacheFilePath(pdbFilePath);
		}
		AtomCache cache = new AtomCache(cacheConfig);

		// Run jobs
		long initialTime = System.nanoTime();
		List<MultipleAlignment> results = new ArrayList<MultipleAlignment>();
		List<SymmetryAxes> axes = new ArrayList<SymmetryAxes>();
		
		try {
			ExecutorService executor = Executors.newFixedThreadPool(threads);
			List<CeSymm> workers = new ArrayList<CeSymm>();
			List<Future<MultipleAlignment>> msaFuture = 
					new ArrayList<Future<MultipleAlignment>>();

			//Run the jobs in parallel Threads
			for(String name: names) {

				Atom[] atoms = null;
				
				try {
					atoms = cache.getAtoms(name);
				} catch (IOException e){
					logger.warn("Could not load structure "+name,e);
					continue;
				} catch (StructureException e) {
					logger.warn("Could not load structure "+name,e);
					continue;
				}

				CeSymm ceSymm = new CeSymm();
				ceSymm.setParameters(params);

				Callable<MultipleAlignment> worker = 
						new CeSymmCallable(atoms, ceSymm);

				workers.add(ceSymm);
				Future<MultipleAlignment> submit = executor.submit(worker);
				msaFuture.add(submit);
			}
			
			for (Future<MultipleAlignment> msa : msaFuture){
				results.add(msa.get());
			}
			for (CeSymm ce : workers){
				axes.add(ce.getSymmetryAxes());
			}
			
			executor.shutdown();

		} catch (InterruptedException e) {
			logger.warn("Job failed: "+e.getMessage());
		} catch (ExecutionException e) {
			logger.warn("Job failed: "+e.getMessage());
		}
		
		long totalSecondsTaken = (System.nanoTime() - initialTime) / 1000000;
		long meanSecondsTaken = (long) (totalSecondsTaken / (float) names.size() / 1000000.0f);
		logger.info("Total runtime: "+totalSecondsTaken + ", mean runtime: "+meanSecondsTaken);

		for (int i=0; i<results.size(); i++){
			
			// Display alignment
			if( displayAlignment ) {
				try {
					SymmetryDisplay.display(results.get(i), axes.get(i));
				} catch (StructureException e) {
					logger.warn("Error displaying symmetry for entry number "+i,e);
				}
			}
	
			// Output alignments
			for(CeSymmWriter writer: writers) {
				try {
					writer.writeAlignment(results.get(i));
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		// close last, in case any writers share output streams
		for(CeSymmWriter writer: writers) {
			writer.close();
		}
	}









	/**
	 * Prompts the user for an input structure using a dialog box
	 * @return The input string, or null if the user cancelled
	 */
	private static String promptUserForStructure() {
		String name;
		// default name
		name = "d1ijqa1";
		//		name = "1G6S";
		name = "1MER.A";
		//		name = "1MER";
		//		name = "1TIM.A";
		//		name = "d1h70a_";
		//name = "2YMS";
		//name = "1HIV";

		name = (String)JOptionPane.showInputDialog(
				null,
				"Structure ID (PDB, SCOP, etc):",
				"Input Structure",
				JOptionPane.PLAIN_MESSAGE,
				null,
				null,
				name);
		if( name != null ) {
			name = name.trim();
		}
		return name;
	}

	/**
	 * Creates the options
	 * @param optionOrder An empty map, which will be filled in with the option
	 *  order
	 * @return all Options
	 */
	@SuppressWarnings("static-access")
	private static Options getOptions(Map<String,Integer> optionOrder) {

		OptionGroup grp;
		Option opt;
		// Note: When adding an option, also add its long name to the optionOrder map
		int optionNum = 0;

		Options options = new Options();
		options.addOption("h","help",false,"Print usage information");
		optionOrder.put("help", optionNum++);
		options.addOption(OptionBuilder.withLongOpt("version")
				.hasArg(false)
				.withDescription("Print CE-Symm version")
				.create());
		optionOrder.put("version", optionNum++);

		// Input file
		options.addOption( OptionBuilder.withLongOpt("input")
				.hasArg(true)
				.withArgName("file")
				.withDescription("File listing whitespace-delimited query structures")
				.create("i"));
		optionOrder.put("input", optionNum++);
		// Output formats
		options.addOption( OptionBuilder.withLongOpt("xml")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output alignment as XML (use --xml=- for standard out).")
				.create("o"));
		optionOrder.put("xml", optionNum++);
		options.addOption( OptionBuilder.withLongOpt("stats")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output a tsv file with the main symmetry info.")
				.create("o"));
		optionOrder.put("stats", optionNum++);
		options.addOption( OptionBuilder.withLongOpt("fasta")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output alignment as FASTA alignment output")
				.create());
		optionOrder.put("fasta", optionNum++);
		options.addOption( OptionBuilder.withLongOpt("fatcat")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output alignment as FATCAT output")
				.create());
		optionOrder.put("fatcat", optionNum++);

		// jmol
		grp = new OptionGroup();
		opt = OptionBuilder
				.withLongOpt("show3d")
				.hasArg(false)
				.withDescription("Force jMol display for each structure "
						+ "[default for <10 structures when specified on command"
						+ " line]")
						.create('j');
		grp.addOption(opt);
		optionOrder.put(opt.getLongOpt(), optionNum++);
		opt = OptionBuilder
				.withLongOpt("noshow3d")
				.hasArg(false)
				.withDescription("Disable jMol display [default with --input "
						+ "or for >=10 structures]")
						.create('J');
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		options.addOptionGroup(grp);

		options.addOptionGroup(grp);
		options.addOption( OptionBuilder.withLongOpt("ordermethod")
				.hasArg(true)
				.withArgName("class")
				.withDescription("Order detection method. Can be a "
						+ "full class name or a short class name from the "
						+ "org.biojava.nbio.structure.align.symmetry.internal package. "
						+ "[default SequenceFunctionOrderDetector]")
						.create());
		optionOrder.put("ordermethod", optionNum++);
		
		options.addOptionGroup(grp);
		options.addOption( OptionBuilder.withLongOpt("refinemethod")
				.hasArg(true)
				.withArgName("class")
				.withDescription("Refiner method. Can be a "
						+ "full class name or a short class name from the "
						+ "org.biojava.nbio.structure.align.symmetry.internal package. "
						+ "[default Single]")
						.create());
		optionOrder.put("refinemethod", optionNum++);
		
		options.addOptionGroup(grp);
		options.addOption( OptionBuilder.withLongOpt("symmtype")
				.hasArg(true)
				.withArgName("class")
				.withDescription("Symmetry Type. Can be a "
						+ "full class name or a short class name from the "
						+ "org.biojava.nbio.structure.align.symmetry.internal package. "
						+ "[default Auto]")
						.create());
		optionOrder.put("symmtype", optionNum++);

		// PDB_DIR
		options.addOption( OptionBuilder.withLongOpt("pdbfilepath")
				.hasArg(true)
				.withArgName("dir")
				.withDescription("Download directory for new "
						+ "structures. Equivalent to passing -DPDB_DIR=dir to the VM. "
						+ "[default temp folder]")
						.create());
		optionOrder.put("pdbfilepath", optionNum++);
		grp = new OptionGroup();
		opt = OptionBuilder
				.withLongOpt("pdbdirsplit")
				.hasArg(false)
				.withDescription("Ignored. For backwards compatibility only. [default]")
				.create();
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		//grp.setSelected(opt);
		opt = OptionBuilder
				.withLongOpt("nopdbdirsplit")
				.hasArg(false)
				.withDescription("Ignored. For backwards compatibility only.")
				.create();
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		options.addOptionGroup(grp);

		options.addOption( OptionBuilder.withLongOpt("threads")
				.hasArg(true)
				.withDescription("Number of threads [default cores-1]")
				.create());
		optionOrder.put("threads", optionNum++);

		//Parameters
		options.addOption( OptionBuilder.withLongOpt("maxgapsize")
				.hasArg(true)
				.withDescription("This parameter configures the maximum gap size "
						+ "G, that is applied during the AFP extension. The "
						+ "larger the value, the longer the calculation time "
						+ "can become, Default value is 30. Set to 0 for no limit.")
						.create()
				);
		optionOrder.put("maxgapsize", optionNum++);

		options.addOption( OptionBuilder.withLongOpt("scoringstrategy")
				.hasArg(true)
				.withDescription("Which scoring function to use: ")
				.create()
				);
		optionOrder.put("scoringstrategy", optionNum++);

		options.addOption( OptionBuilder.withLongOpt("winsize")
				.hasArg(true)
				.withDescription("This configures the fragment size m of Aligned Fragment Pairs (AFPs).")
				.create()
				);
		optionOrder.put("winsize", optionNum++);

		options.addOption( OptionBuilder.withLongOpt("maxrmsd")
				.hasArg(true)
				.withDescription("The maximum RMSD at which to stop alignment "
						+ "optimization. (default: unlimited=99)")
						.create()
				);
		optionOrder.put("maxrmsd", optionNum++);

		options.addOption( OptionBuilder.withLongOpt("gapopen")
				.hasArg(true)
				.withDescription("Gap opening penalty during alignment optimization [default: 5.0].")
				.create()
				);
		optionOrder.put("gapopen", optionNum++);

		options.addOption( OptionBuilder.withLongOpt("gapextension")
				.hasArg(true)
				.withDescription("Gap extension penalty during alignment optimization [default: 0.5].\n")
				.create()
				);
		optionOrder.put("gapextension", optionNum++);
		
		options.addOption( OptionBuilder.withLongOpt("multaxes")
				.hasArg(true)
				.withDescription("Run iteratively the algorithm to find multiple symmetry levels [default: true].\n")
				.create()
				);
		optionOrder.put("multaxes", optionNum++);
		
		options.addOption( OptionBuilder.withLongOpt("opt")
				.hasArg(true)
				.withDescription("Optimize the resulting symmetry alignment [default: true].\n")
				.create()
				);
		optionOrder.put("opt", optionNum++);
		
		options.addOption( OptionBuilder.withLongOpt("threshold")
				.hasArg(true)
				.withDescription("The symmetry threshold. TM-scores above this value"
						+ "will be considered significant results [default: 0.5, interval [0.0,1.0]].\n")
				.create()
				);
		optionOrder.put("threshold", optionNum++);
		
		options.addOption( OptionBuilder.withLongOpt("maxorder")
				.hasArg(true)
				.withDescription("The maximum number of symmetric subunits [default: 8].\n")
				.create()
				);
		optionOrder.put("maxorder", optionNum++);
		
		options.addOption( OptionBuilder.withLongOpt("rndseed")
				.hasArg(true)
				.withDescription("The random seed used in optimization, for reproducibility"
						+ "of the results [default: 0].\n")
				.create()
				);
		optionOrder.put("rndseed", optionNum++);
		
		options.addOption( OptionBuilder.withLongOpt("minlen")
				.hasArg(true)
				.withDescription("The minimum length, expressed in number of core "
						+ "aligned residues, of a symmetric subunit [default: 15].\n")
				.create()
				);
		optionOrder.put("minlen", optionNum++);
		
		options.addOption( OptionBuilder.withLongOpt("dcutoff")
				.hasArg(true)
				.withDescription("The maximum distance, in A, allowed between any two aligned "
						+ "residue positions [default: 7.0].\n")
				.create()
				);
		optionOrder.put("dcutoff", optionNum++);

		options.addOption( OptionBuilder.withLongOpt("scopversion")
				.hasArg(true)
				.withArgName("version")
				.withDescription("Version of SCOP or SCOPe to use "
						+ "when resolving SCOP identifiers [defaults to latest SCOPe]")
						.create());
		optionOrder.put("scopversion", optionNum++);

		return options;
	}

	/**
	 * Parse a whitespace-delimited file containing structure names
	 * @throws FileNotFoundException 
	 */
	public static List<String> parseInputStructures(String filename) throws FileNotFoundException {
		File file = new File(filename);
		Scanner s = new Scanner(file);

		List<String> structures = new ArrayList<String>();
		while(s.hasNext()) {
			String name = s.next();
			if(name.startsWith("#")) {
				//comment
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
	 * Parent class for all output formats
	 * All methods are empty stubs, which should be overridden to write
	 * data to the writer.
	 */
	private static abstract class CeSymmWriter {
		protected PrintWriter writer;
		public CeSymmWriter(PrintWriter writer) {
			this.writer = writer;
		}
		public CeSymmWriter(String filename) throws IOException {
			this(openOutputFile(filename));
		}
		public void writeHeader() throws IOException {}
		public void writeAlignment(MultipleAlignment alignment) throws IOException {}
		public void close() {
			if(writer != null) {
				writer.flush();
				writer.close();
			}
		}
		
		/**
		 * Opens 'filename' for writing.
		 * @param filename Name of output file, or '-' for standard out
		 * @throws IOException 
		 */
		public static PrintWriter openOutputFile(String filename) throws IOException {
			if(filename.equals("-")) {
				return new PrintWriter(System.out,true);
			}
			return new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		}
	}
	
	private static class XMLWriter extends CeSymmWriter {
		public XMLWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeAlignment(MultipleAlignment msa) throws IOException {
			writer.append(MultipleAlignmentWriter.toXML(msa.getEnsemble()));
			writer.flush();
		}
	}
	
	private static class FatcatWriter extends CeSymmWriter {
		public FatcatWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeAlignment(MultipleAlignment msa) {
			writer.write(MultipleAlignmentWriter.toFatCat(msa));
			writer.println("//");
			writer.flush();
		}
	}
	private static class FastaWriter extends CeSymmWriter {
		public FastaWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeAlignment(MultipleAlignment msa) {
			writer.write(MultipleAlignmentWriter.toFASTA(msa));
			writer.println("//");
			writer.flush();
		}
	}
	private static class StatsWriter extends CeSymmWriter {
		
		public StatsWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeHeader() {
			writer.println("Name\t" +
					"Order\t" +
					"PG\t" +
					"isRefined\t" +
					"Avg-TMscore\t" +
					"Avg-RMSD\t" +
					"SubunitLength\t" +
					"Length\t" +
					"CoreLength\t" +
					"Coverage\t");
			writer.flush();
		}
		@Override
		public void writeAlignment(MultipleAlignment msa)
				throws IOException {
			
			String pg = null;
			int order = 1;
			MultipleAlignment full = msa;
			int subunitLen = 0;
			boolean refined = SymmetryTools.isRefined(msa);
			
			if (refined){
				full = SymmetryTools.toFullAlignment(msa);
				pg = SymmetryTools.getQuaternarySymmetry(msa).getSymmetry();
				order = msa.size();
				subunitLen = msa.length();
			}
			double structureLen = msa.getEnsemble().getAtomArrays().get(0).length;
			double coverage = full.length() / structureLen;
			
			writer.format("%s\t%d\t%s\t%b\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f\n",
					msa.getEnsemble().getStructureNames().get(0),
					order,
					pg,
					refined,
					msa.getScore(MultipleAlignmentScorer.AVGTM_SCORE),
					msa.getScore(MultipleAlignmentScorer.RMSD),
					subunitLen,
					full.length(),
					full.getCoreLength(),
					coverage
					);
			writer.flush();
		}
	}

	public static class CeSymmCallable implements Callable<MultipleAlignment> {

		private Atom[] atoms;
		private CeSymm ceSymm;

		public CeSymmCallable(Atom[] atoms, CeSymm ceSymm){
			this.atoms = atoms;
			this.ceSymm = ceSymm;
		}

		@Override
		public MultipleAlignment call() throws Exception {
			return ceSymm.analyze(atoms);
		}
	}

}
