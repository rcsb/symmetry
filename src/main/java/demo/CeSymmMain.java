package demo;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.order.OrderDetectionFailedException;
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
		PrintStream xmlOut = null;
		PrintStream fatcatOut = null;
		PrintStream ceOut = null;
		PrintStream tsvOut = null;
		if(cli.hasOption("xml")) {
			xmlOut = openOutputFile(cli.getOptionValue("xml"));
		}
		if(cli.hasOption("fatcat")) {
			fatcatOut = openOutputFile(cli.getOptionValue("fatcat"));
		}
		if(cli.hasOption("ce")) {
			ceOut = openOutputFile(cli.getOptionValue("ce"));
		}
		if(cli.hasOption("tsv")) {
			tsvOut = openOutputFile(cli.getOptionValue("tsv"));
		} else if(cli.hasOption("verbose")) {
			tsvOut = openOutputFile("-");
		}
		// TODO pdbOut (needs individual files)

		// Done parsing arguments

		UserConfiguration cacheConfig = new UserConfiguration();
		if(pdbFilePath != null && !pdbFilePath.isEmpty()) {
			cacheConfig.setPdbFilePath(pdbFilePath);
			cacheConfig.setCacheFilePath(pdbFilePath);
		}
		if(pdbDirSplit != null) {
			cacheConfig.setSplit(pdbDirSplit);
		}
		AtomCache cache = new AtomCache(cacheConfig);

		for(String name: names) {
			try {

				// Perform alignment to determine axis
				Atom[] ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name,null,cache));
				Atom[] ca2 = StructureTools.cloneCAArray(ca1);
				CeSymm ce = new CeSymm();
				AFPChain alignment = ce.align(ca1, ca2);
				alignment.setName1(name);
				alignment.setName2(name);
				RotationAxis axis = new RotationAxis(alignment);

				// Display alignment
				if( displayAlignment ) {
					StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(alignment, ca1, ca2);
					jmol.evalString(axis.getJmolScript(ca1));
				}

				// Significance
				boolean significant = false;
				if(useOrder) {
					significant = ce.isSignificant();
				} else {
					//TODO remove hard coded threshold
					significant = alignment.getTMScore() >= .4;
				}

				// Order
				int symmNr = 1;
				if( significant ) {
					try {
						symmNr = detector.calculateOrder(alignment, ca1);
					} catch (OrderDetectionFailedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}

				// Print result
				System.out.format("%s\tTMscore %f\tOrder %d\tSignificant %s%n",name,alignment.getTMScore(),symmNr,significant?"Y":"N");
				// Output
				// TODO implement output
				System.out.println(alignment.toFatcat(ca1,ca2));

			} catch (IOException e) {
				e.printStackTrace();
			} catch (StructureException e) {
				e.printStackTrace();
			}
		}

		// Close outpus
		if(xmlOut != null) {
			xmlOut.close();
		}
		if(fatcatOut != null) {
			fatcatOut.close();
		}
		if(ceOut != null) {
			ceOut.close();
		}
		if(tsvOut != null) {
			tsvOut.close();
		}
	}

	/**
	 * Opens 'filename' for writing.
	 * @param filename Name of output file, or '-' for standard out
	 */
	private static PrintStream openOutputFile(String filename) {
		if(filename.equals("-")) {
			return System.out;
		}
		// TODO Auto-generated method stub
		return System.out;
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
		name = "1HIV";

		name = (String)JOptionPane.showInputDialog(
				null,
				"Structure ID (PDB, SCOP, etc):",
				"Input Structure",
				JOptionPane.PLAIN_MESSAGE,
				null,
				null,
				name);
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
				.withDescription("Output alignment as XML (use --xml=- for standard out")
				.create("o"));
		optionOrder.put("xml", optionNum++);
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
		//TODO pdb output
		//		options.addOption( OptionBuilder.withLongOpt("pdb")
		//				.hasArg(true)
		//		.withArgName("file")
		//				.withDescription("Output alignment as two-model PDB file")
		//				.create());
		options.addOption( OptionBuilder.withLongOpt("tsv")
				.hasArg(true)
				.withArgName("file")
				.withDescription("Output alignment as tab-separated file")
				.create());
		optionOrder.put("tsv", optionNum++);


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
}
