package demo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.IllegalFormatException;
import java.util.List;
import java.util.Map;
import java.util.MissingFormatArgumentException;
import java.util.Scanner;

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
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.gui.DisplayAFP;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.io.util.FileDownloadUtils;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.census3.AdditionalScoreList;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusResultList;
import org.biojava3.structure.align.symm.census3.CensusScoreList;
import org.biojava3.structure.align.symm.census3.MapScoreList;
import org.biojava3.structure.align.symm.census3.run.Census;
import org.biojava3.structure.align.symm.census3.run.CensusJob;
import org.biojava3.structure.align.symm.order.OrderDetector;
import org.biojava3.structure.align.symm.order.SequenceFunctionOrderDetector;

/**
 * Main executable for running CE-Symm
 * 
 * Run with -h for usage, or without arguments for interactive mode
 * @author spencer
 *
 */
public class CeSymmMain {


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

		CommandLine cli;
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

		// Use order? [default true]
		boolean useOrder = cli.hasOption("order") || !cli.hasOption("noorder");

		// Order method
		String orderMethod = null;
		if(cli.hasOption("ordermethod") ) {
			orderMethod = cli.getOptionValue("ordermethod");
		}
		OrderDetector detector = createOrderDetector(orderMethod);
		if(detector == null) {
			System.exit(1);
		}

		// AtomCache options
		String pdbFilePath = null;
		if( cli.hasOption("pdbfilepath") ) {
			pdbFilePath = cli.getOptionValue("pdbfilepath");
			pdbFilePath = FileDownloadUtils.expandUserHome(pdbFilePath);
		}
		Boolean pdbDirSplit = null;
		if( cli.hasOption("nopdbdirsplit") ) {
			pdbDirSplit = false;
		}
		if( cli.hasOption("pdbdirsplit") ) {
			pdbDirSplit = true;
		}

		//TODO threads
		//		Integer threads = null;
		//		if(cli.hasOption("threads")) {
		//			threads = new Integer(cli.getOptionValue("threads"));
		//		}

		// SCOP version
		if( cli.hasOption("scopversion")) {
			String scopVersion = cli.getOptionValue("scopversion");

			//TODO add validation to version 
			ScopFactory.setScopDatabase(scopVersion);
		}

		// Output formats
		List<CeSymmWriter> writers = new ArrayList<CeSymmWriter>();

		if(cli.hasOption("stats")) {
			String filename = cli.getOptionValue("stats");
			try {
				writers.add(new SimpleWriter(filename,useOrder?detector:null));
			} catch (IOException e) {
				System.err.println("Error: Ignoring file "+filename+".");
				System.err.println(e.getMessage());
			}
		}
		else if(!cli.hasOption("nostats")) {
			try {
				//default stdout output
				writers.add(new SimpleWriter("-",useOrder?detector:null));
			} catch (IOException e1) {}
		}

		if(cli.hasOption("xml")) {
			String filename = cli.getOptionValue("xml");
			try {
				writers.add(new XMLWriter(filename));
				//writers.add(new AltXMLWriter(filename));
			} catch (IOException e) {
				System.err.println("Error: Ignoring file "+filename+".");
				System.err.println(e.getMessage());
			}
		}
		if(cli.hasOption("html")) {
			String filename = cli.getOptionValue("html");
			try {
				writers.add(new HTMLWriter(filename));
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
		if(cli.hasOption("ce")) {
			String filename = cli.getOptionValue("ce");
			try {
				writers.add(new CeWriter(filename));
			} catch (IOException e) {
				System.err.println("Error: Ignoring file "+filename+".");
				System.err.println(e.getMessage());
			}
		}

		if(cli.hasOption("tsv")) {
			String filename = cli.getOptionValue("tsv");
			try {
				writers.add(new TSVWriter(filename));
			} catch (IOException e) {
				System.err.println("Error: Ignoring file "+filename+".");
				System.err.println(e.getMessage());
			}
		} else if(cli.hasOption("verbose")) {
			try {
				writers.add(new TSVWriter("-"));
			} catch (IOException e) {} //shouldn't happen
		}

		if(cli.hasOption("pdb")) {
			String filespec = cli.getOptionValue("pdb");
			writers.add(new PDBWriter(filespec));
		}

		// Done parsing arguments

		// Configure atomcache
		UserConfiguration cacheConfig = new UserConfiguration();
		if(pdbFilePath != null && !pdbFilePath.isEmpty()) {
			cacheConfig.setPdbFilePath(pdbFilePath);
			cacheConfig.setCacheFilePath(pdbFilePath);
		}
		if(pdbDirSplit != null) {
			cacheConfig.setSplit(pdbDirSplit);
		}
		AtomCache cache = new AtomCache(cacheConfig);

		// Add as option to background alignment GUI
		StructureAlignmentFactory.addAlgorithm(new CeSymm());

		CensusResultList results = new CensusResultList();

		//print headers
		for(CeSymmWriter writer: writers) {
			try {
				writer.writeHeader();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		// Run jobs
		long totalTimeTaken = 0;
		for(String name: names) {
			try {

				final CensusJob calc = CensusJob.setUpJob(name, 1, Census.AlgorithmGiver.getDefault(),
						Census.getDefaultAfpChainCensusRestrictor(), cache);
				calc.setRecordAlignmentMapping(true);
				calc.setStoreAfpChain(true);
				calc.setOrderDetector(detector);

				CensusResult result = calc.call();
				if (result == null) continue; 
				totalTimeTaken += calc.getTimeTaken();

				results.add(result);
				// Probably an abuse of this property, but I'm not sure how it
				// was intended. Used by SimpleWriter
				Map<String, Number> scoreMap = new HashMap<String, Number>();
				scoreMap.put( "timeMillis",calc.getTimeTaken());
				AdditionalScoreList moreScores = new MapScoreList( scoreMap );
				result.getScoreList().setAdditionalScoreList(moreScores);

				results.setMeanSecondsTaken((float) (totalTimeTaken / (float) names.size() / 1000.0f));

				// Perform alignment to determine axis
				Atom[] ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(result.getId(),null,cache));
				Atom[] ca2 = StructureTools.cloneCAArray(ca1);
				AFPChain alignment = calc.getAfpChain();
				//alignment.setName1(name);
				//alignment.setName2(name);
				RotationAxis axis = result.getAxis().toRotationAxis();

				// Display alignment
				if( displayAlignment ) {
					StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(alignment, ca1, ca2);
					jmol.evalString(axis.getJmolScript(ca1));
				}

				// Output alignments
				for(CeSymmWriter writer: writers) {
					try {
						writer.writeAlignment(alignment, result, ca1, ca2);
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			} catch (IOException e) {
				e.printStackTrace();
			} catch (StructureException e) {
				e.printStackTrace();
			}
		}

		// Output footers
		for(CeSymmWriter writer: writers) {
			try {
				writer.writeFooter(results);
			} catch (IOException e) {
				e.printStackTrace();
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
		options.addOption( OptionBuilder.withLongOpt("html")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output alignment as HTML output")
				.create());
		optionOrder.put("html", optionNum++);
		options.addOption( OptionBuilder.withLongOpt("ce")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output alignment as CE output")
				.create());
		optionOrder.put("ce", optionNum++);
		options.addOption( OptionBuilder.withLongOpt("fatcat")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output alignment as FATCAT output")
				.create());
		optionOrder.put("fatcat", optionNum++);
		options.addOption( OptionBuilder.withLongOpt("pdb")
				.hasArg(true)
				.withArgName("dir")
				.withDescription("Output each alignment as a two-model PDB file. "
						+ "The argument may be a directory or a formatting string, "
						+ "where \"%s\" will be replaced with the structure name. "
						+ "[default \"%s.cesymm.pdb\"]")
						.create());
		optionOrder.put("pdb", optionNum++);
		options.addOption( OptionBuilder.withLongOpt("tsv")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output alignment as tab-separated file")
				.create());
		optionOrder.put("tsv", optionNum++);
		options.addOption( OptionBuilder.withLongOpt("stats")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output tab-separated file giving alignment "
						+ "statistics, one line per structure [defaults to stdout,"
						+ " unless -q is given]")
				.create());
		optionOrder.put("stats", optionNum++);



		options.addOption(OptionBuilder.withLongOpt("nostats")
				.hasArg(false)
				.withDescription("Do not output default statistics to standard "
						+ "out  (equivalent to \"--stats=/dev/null\")")
				.create('q'));
		optionOrder.put("nostats", optionNum++);
		options.addOption(OptionBuilder.withLongOpt("verbose")
				.hasArg(false)
				.withDescription("Print detailed output (equivalent to \"--tsv=-\")")
				.create('v'));
		optionOrder.put("verbose", optionNum++);

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

		// order
		grp = new OptionGroup();
		opt = OptionBuilder
				.withLongOpt("order")
				.hasArg(false)
				.withDescription("Use TM-Score with order for deciding significance. [default]")
				.create('t');
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		//grp.setSelected(opt);
		opt = OptionBuilder
				.withLongOpt("noorder")
				.hasArg(false)
				.withDescription("Use TM-Score alone for deciding significance.")
				.create('T');
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		options.addOptionGroup(grp);
		options.addOption( OptionBuilder.withLongOpt("ordermethod")
				.hasArg(true)
				.withArgName("class")
				.withDescription("Order detection method. Can be a "
						+ "full class name or a short class name from the "
						+ "org.biojava3.structure.align.symm.order package. "
						+ "[default SequenceFunctionOrderDetector]")
						.create());
		optionOrder.put("ordermethod", optionNum++);

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
				.withDescription("Indicates that --pdbfilepath is split into "
						+ "multiple subdirs, like the ftp site. [default]")
						.create();
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		//grp.setSelected(opt);
		opt = OptionBuilder
				.withLongOpt("nopdbdirsplit")
				.hasArg(false)
				.withDescription("Indicates that --pdbfilepath should be a single directory.")
				.create();
		optionOrder.put(opt.getLongOpt(), optionNum++);
		grp.addOption(opt);
		options.addOptionGroup(grp);

		// misc
		//TODO threads
		//		options.addOption( OptionBuilder.withLongOpt("threads")
		//				.hasArg(true)
		//				.withDescription("Number of threads [default cores-1]")
		//				.create());
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


	/**
	 * Creates an OrderDetector from a class name.
	 * 
	 * <p>
	 * Accepts the following inputs:<ol>
	 *  <li>null or "" returns a SequenceFunctionOrderDetector
	 *  <li>The full class path to a class implementing OrderDetector and containing a default constructor
	 *  <li>The class name for any class in the org.biojava3.structure.align.symm.order package.
	 * </ol>
	 * 
	 * @param method Name of the OrderDetector method
	 * @return An OrderDetector instance, or null for invalid input
	 */
	private static OrderDetector createOrderDetector(String method) {	

		if(method == null || method.isEmpty()) {
			return new SequenceFunctionOrderDetector();
		}

		ClassLoader cl = CeSymmMain.class.getClassLoader();
		Class<?> klass = null;
		// try full class name
		try {
			klass = cl.loadClass(method);
		} catch( ClassNotFoundException e) {
			// ignore
		}

		// try order package
		try {
			String fullname = OrderDetector.class.getPackage().getName()+"."+method;
			klass = cl.loadClass(fullname);
		} catch( ClassNotFoundException e) {
			//ignore
		}

		// Give up if that didn't work
		if(klass == null) {
			System.err.format("Error: Method '%s' not found.%n",method);
			return null;
		}

		// Instantiate default constructor
		OrderDetector detector = null;
		try {
			Constructor<?> constructor = klass.getConstructor();
			detector = (OrderDetector) constructor.newInstance();
		} catch (ClassCastException e) {
			// Not an OrderDetector
			System.err.println("Error: "+method+" is not an OrderDetector.");
		} catch( NoSuchMethodException e) {
			// No default constructor
			System.err.println("Error: Unable to use "+method+" because it lacks a default constructor");
		} catch (IllegalArgumentException e) {
			// Shouldn't happen–bad argument types
			System.err.println("Error: [Bug] Error with constructor arguments to "+method);
		} catch (InstantiationException e) {
			// Abstract class
			System.err.println("Error: Can't instantiate abstract class "+method);
		} catch (IllegalAccessException e) {
			// constructor is private
			System.err.println("Error: "+method+" lacks a public default constructor");
		} catch (InvocationTargetException e) {
			// Constructor threw an exception
			System.err.println("Error: Exception while creating "+method);
			e.getCause().printStackTrace();
		}

		return detector;
	}

	// Output formats

	/**
	 * Parent class for all output formats
	 * All methods are empty stubs, which should be overridden to write
	 * data to the writer.
	 * @author blivens
	 *
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
		public void writeAlignment(AFPChain alignment, CensusResult result, Atom[] ca1, Atom[] ca2) throws IOException {}
		public void writeFooter(CensusResultList results) throws IOException {}
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
		public void writeFooter(CensusResultList results) throws IOException {
			writer.write(results.toXML());
			writer.flush();
		}
	}
	private static class AltXMLWriter extends CeSymmWriter {
		public AltXMLWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeAlignment(AFPChain afpChain, CensusResult results,
				Atom[] ca1, Atom[] ca2) throws IOException {
			writer.append(AFPChainXMLConverter.toXML(afpChain, ca1, ca2));
			writer.flush();
		}
	}
	private static class HTMLWriter extends CeSymmWriter {
		public HTMLWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeHeader() {
			StringBuilder html = new StringBuilder();
			html.append("<html>\n");
			html.append("  <head>\n");
			html.append("    <title>CeSymm</title>\n");
			html.append("    <style type=\"text/css\">\n" +
					".aligmentDisp span { color:#999999; }\n" + 
					"span.m { color:#008800; }\n" + 
					"span.dm, #alignment span.na { color: #A66A00; } \n" + 
					"span.sm { color: #D460CF; }\n" + 
					"span.qg, #alignment span.hg { color: #104BA9; }\n" + 
					".aligmentDisp { border-style: solid; border-width: 1px; font-family: \"Courier New\", Courier, monospace; font-size: 1.1em; margin: 10px; overflow-x:auto; padding:10px; width:90%; white-space: nowrap; }\n" + 
					"\n" + 
					".alignmentBlock11 { background-color: rgb(255, 165, 0); }\n" + 
					".alignmentBlock12 { background-color: rgb(255, 229, 0); }\n" + 
					".alignmentBlock13 { background-color: rgb(217, 255, 0); }\n" + 
					".alignmentBlock14 { background-color: rgb(154, 255, 0); }\n" + 
					".alignmentBlock15 { background-color: rgb( 90, 255, 0); }\n" + 
					"\n" + 
					".alignmentBox11   { background-color: rgb(255, 165, 0);  border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					".alignmentBox12   { background-color: rgb(255, 229, 0);  border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					".alignmentBox13   { background-color: rgb(217, 255, 0);  border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					".alignmentBox14   { background-color: rgb(154, 255, 0);  border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					".alignmentBox15   { background-color: rgb( 90, 255, 0);  border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					"\n" + 
					".alignmentBlock21 { background-color: rgb(0, 255, 255); }\n" + 
					".alignmentBlock22 { background-color: rgb(0, 178, 255); }\n" + 
					".alignmentBlock23 { background-color: rgb(0, 102, 255); }\n" + 
					".alignmentBlock24 { background-color: rgb(0,  25, 255); }\n" + 
					".alignmentBlock25 { background-color: rgb(51, 0, 255);  } \n" + 
					"\n" + 
					".alignmentBox21   { background-color: rgb(0, 255, 255); border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					".alignmentBox22   { background-color: rgb(0, 178, 255); border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					".alignmentBox23   { background-color: rgb(0, 102, 255); border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					".alignmentBox24   { background-color: rgb(0,  25, 255); border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					".alignmentBox25   { background-color: rgb(51,  0, 255); border-color: black; border-style: solid; border-width: 1px; margin: 5px 10px 10px 0pt; height: 10px; width: 10px; }\n" + 
					"    </style>\n");
			html.append("  </head>\n");
			html.append("  <body>\n");
			html.append("    <div id=\"mainContent\">\n");
			writer.append(html.toString());
			writer.flush();
		}
		@Override
		public void writeAlignment(AFPChain alignment, CensusResult result, Atom[] ca1, Atom[] ca2) {
			writer.write("<h3>"+alignment.getName1()+"</h3>\n");
			writer.write("<pre>\n");
			writer.write(AfpChainWriter.toWebSiteDisplay(alignment, ca1, ca2));
			writer.write("</pre>\n");
			writer.flush();
		}
		@Override
		public void writeFooter(CensusResultList results) {
			StringBuilder html = new StringBuilder();
			html.append("    </div>");
			html.append("  </body>\n");
			html.append("</html>\n");
			writer.append(html.toString());
			writer.flush();
		}

	}
	private static class FatcatWriter extends CeSymmWriter {
		public FatcatWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeAlignment(AFPChain alignment, CensusResult result, Atom[] ca1, Atom[] ca2) {
			writer.write(alignment.toFatcat(ca1,ca2));
			writer.println("//");
			writer.flush();
		}
	}
	private static class CeWriter extends CeSymmWriter {
		public CeWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeAlignment(AFPChain alignment, CensusResult result, Atom[] ca1, Atom[] ca2) {
			writer.write(alignment.toCE(ca1,ca2));
			writer.println("//");
			writer.flush();
		}
	}
	private static class TSVWriter extends CeSymmWriter {
		public TSVWriter(String filename) throws IOException {
			super(filename);
		}
		@Override
		public void writeAlignment(AFPChain alignment, CensusResult result, Atom[] ca1, Atom[] ca2) {
			writer.write(AfpChainWriter.toAlignedPairs(alignment, ca1, ca2));
			writer.println("//");
			writer.flush();
		}
	}
	private static class PDBWriter extends CeSymmWriter {
		private String pdbFormat;
		public PDBWriter(String format) {
			super((PrintWriter)null);
			pdbFormat = createPDBFilenameFormatString(format);
		}
		/**
		 * @param pdbFormat
		 * @param filespec
		 * @return
		 */
		private static String createPDBFilenameFormatString( String filespec) {
			String pdbFormat;
			if(filespec.isEmpty()) {
				filespec = ".";
			}
			File pdbOutDir = new File(filespec);
			// 1) It is a directory. Use default filenames
			if(pdbOutDir.isDirectory()) { //should also cover ""
				pdbFormat = new File(pdbOutDir,"%s.cesymm.pdb").getPath();
			} else {
				// Split into file/directory & check for existence
				String name = pdbOutDir.getName();
				assert(!name.isEmpty());

				String parent = pdbOutDir.getParent();
				if(parent == null) parent = ".";
				pdbOutDir = new File(parent);
				if(!pdbOutDir.isDirectory() || !pdbOutDir.canWrite()) {
					System.err.println("Error: unable to write to "+filespec);
					System.exit(1);
					return null;
				}
				try {
					String.format(filespec); // throws exception if filespec has "%s"
					// 2) Not a format string. Append structure name and extension
					pdbFormat = filespec+".%s.pdb";
				} catch(MissingFormatArgumentException e) {
					// 3) The filename contains %s, making it a valid format string
					try {
						String.format(filespec,"");
						pdbFormat = filespec;
					} catch(IllegalFormatException f) {
						System.err.println("Illegal pdb formating string (use a single %s): "+filespec);
						System.exit(1);
						return null;
					}
				}
			}
			return pdbFormat;
		}
		@Override
		public void writeAlignment(AFPChain alignment, CensusResult result, Atom[] ca1, Atom[] ca2) throws IOException {
			try {
				Structure s = DisplayAFP.createArtificalStructure(alignment, ca1, ca2);
				
				// If the input was from a file then we need to remove the path
				String escapedName = alignment.getName1();
				File file = new File(FileDownloadUtils.expandUserHome(escapedName));
				if(file.exists()) {
					escapedName = file.getName();
				}
				escapedName = escapedName.replaceAll(File.separator, "_");
				String filename = String.format(pdbFormat,escapedName);
				PrintWriter pdbOut = openOutputFile(filename);

				String pdb = s.toPDB();
				pdbOut.write(pdb);
				pdbOut.flush();
				pdbOut.close();
			} catch(StructureException e) {
				e.printStackTrace();
			}
		}
	}
	private static class SimpleWriter extends CeSymmWriter {
		OrderDetector detector;
		/**
		 * Significance calculations are currently hard-coded (TODO) to be
		 * use a hard-coded TM-Score threshold if detector is null, or use
		 * order plus the TM-Score otherwise.
		 * @param detector Either an OrderDetector or null
		 */
		public SimpleWriter(String filename,OrderDetector detector) throws IOException {
			super(filename);
			this.detector = detector;
		}
		@Override
		public void writeHeader() {
			writer.println("Name\t" +
					"Symm\t" +
					"MinOrder\t" +
					"TMscore\t" +
					"ZScore\t" +
					"CEScore\t" +
					"PValue\t" +
					"RMSD\t" +
					"Length\t" +
					"Coverage\t" +
					"%ID\t" +
					"%Sim\t" +
					"time");
			writer.flush();
		}
		@Override
		public void writeAlignment(AFPChain alignment, CensusResult result, Atom[] ca1, Atom[] ca2)
				throws IOException {
			CensusScoreList scores = result.getScoreList();

			//TODO Use a consistent method for significance, rather than hard coding in output(!) code
			boolean significant = false;
			if( detector != null ) {
				try {
					significant = CeSymm.isSignificant(alignment, detector, ca1);
				} catch (StructureException e) {
					e.printStackTrace();
				}
			} else {
				significant = scores.getTmScore() >= 0.4;
			}

			Number timeMillis = scores.getAdditionalScoreList().getScore("timeMillis");
			writer.format("%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.1f\t%.1f\t%.4f%n",
					alignment.getName1(),
					significant?"Y":"N",
					significant?result.getOrder():1,
					scores.getTmScore(),
					scores.getzScore(),
					alignment.getAlignScore(),
					alignment.getProbability(),
					scores.getRmsd(),
					scores.getAlignLength(),
					scores.getIdentity()*100,
					scores.getSimilarity()*100,
					timeMillis.longValue()/1000.
					);
			writer.flush();
		}
	}
}
