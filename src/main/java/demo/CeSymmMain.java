package demo;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import javax.swing.JOptionPane;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.RotationAxis;
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

	/**
	 * Stub order detector. Always returns 100.
	 * 
	 * @author blivens
	 *
	 */
	public class AlwaysSymmetricOrderDetector implements OrderDetector {

		/**
		 * @return 100 for all cases.
		 */
		@Override
		public int calculateOrder(AFPChain afpChain, Atom[] ca)
				throws OrderDetectionFailedException {
			return 100;
		}

	}

	public static void main(String[] args) {
		// Begin argument parsing
		final String usage = "[OPTIONS] [structures...]";
		final String header = "Determine the order for each structure, which may " +
				"be PDB IDs, SCOP domains, or file paths. If none are given, the " +
				"user will be prompted at startup.";
		Options options = getOptions();
		CommandLineParser parser = new GnuParser();
		HelpFormatter help = new HelpFormatter();

		CommandLine cli;
		try {
			cli = parser.parse(options,args,false);
			if(cli.hasOption('h')) {
				help.printHelp(usage, header, options, "");
				System.exit(1);
				return;
			}
		} catch (ParseException e) {
			System.err.println("Error: "+e.getMessage());
			help.printHelp(usage, header, options, "");
			System.exit(1);
			return;
		}

		args = cli.getArgs();


		String[] names;
		if(args.length == 0) {
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

			if( name == null) {
				//cancel
				return;
			}
			names = new String[] {name};
		} else {
			// take names from the command line arguments
			names = args;
		}


		// Show jmol?
		boolean displayAlignment = cli.hasOption('j') || !cli.hasOption('J');
		// Show alignment?
		boolean printAlignment = cli.hasOption('a') && !cli.hasOption('A');
		
		// Symmetry method
		String method = null;
		if(cli.hasOption('m') ) {
			method = cli.getOptionValue('m');
		}
		OrderDetector detector = createOrderDetector(method);
		if(detector == null) {
			System.exit(1);
		}
		
		// Significance Method
		boolean useOrder = cli.hasOption('w') || !cli.hasOption('W');
		// Done parsing arguments


		for(String name: names) {
			try {

				// Perform alignment to determine axis
				Atom[] ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
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
				// Print alignment
				if(printAlignment)
					System.out.println(alignment.toFatcat(ca1,ca2));

			} catch (IOException e) {
				e.printStackTrace();
			} catch (StructureException e) {
				e.printStackTrace();
			}
		}
	}

	private static Options getOptions() {
		Options options = new Options();
		options.addOption("h","help",false,"print help");
		options.addOption("j","jmol",false,"enable jMol display. [default]");
		options.addOption("J","nojmol",false,"disable jMol display.");
		options.addOption("a","alignment",false,"print alignment.");
		options.addOption("A","noalignment",false,"don't print alignment [default]");
		options.addOption("d","detector",true,"Order detection method. Can be a full class name,"
				+ "a short class name from the org.biojava3.structure.align.symm.order package,"
				+ "or one of the shortcuts 'tm' or 'order'. [default 'tm']");
		options.addOption("w","withorder",false,"Use TM-Score with order for deciding "
				+ "signifiance. [default]");
		options.addOption("W","withoutorder",false,"Use TM-Score alone for deciding "
				+ "signifiance.");
		return options;
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
