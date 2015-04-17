package demo;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JOptionPane;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava.nbio.structure.align.symm.order.OrderDetector;
import org.biojava.nbio.structure.align.symm.order.RotationOrderDetector;
import org.biojava.nbio.structure.align.symm.order.RotationOrderDetector.RotationOrderMethod;
import org.biojava.nbio.structure.align.util.RotationAxis;

public class DetectOrder {

	/**
	 * Rotates a structure around an axis by regular intervals and prints the
	 * {@link RotationOrderDetector#superpositionDistance(Atom[], Atom[]) superpositionDistance}
	 * for each.
	 * 
	 * Output is two tab-delimited columns giving the angle and distance.
	 * @param ca Structure to consider
	 * @param axis Axis to rotate around
	 * @param angleIncr Angle, in radians, between samples
	 * @param out Output stream
	 * @throws StructureException If the atoms can't be cloned or manipulated.
	 */
	public static void printSuperpositionDistance(Atom[] ca, RotationAxis axis,double angleIncr, PrintStream out) throws StructureException {
		final int steps = (int)Math.floor(Math.PI/angleIncr);

		Atom[] ca2 = StructureTools.cloneCAArray(ca);

		for (int step=0; step<steps;step++) {
			double dist = RotationOrderDetector.superpositionDistance(ca, ca2);
			double angle = angleIncr*step;

			// Rotate for next step
			axis.rotate(ca2, angleIncr);

			out.format("%f\t%f%n",angle,dist);
		}
	}

	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		// Begin argument parsing
		final String usage = "[OPTIONS] [structure]";
		final String header = "Determine the order for <structure>, which may " +
				"be a PDB ID, SCOP domain, or file path. If none is given, the " +
				"user will be prompted at startup.";

		Options options = new Options();
		options.addOption("h","help",false,"print help");
		options.addOption(OptionBuilder.withArgName("int")
				.hasArg()
				.withLongOpt("max-order")
				.withType(Number.class)
				.withDescription("maximum order to consider [default 8]")
				.create("x") );
		options.addOption(OptionBuilder.withArgName("degrees")
				.hasArg()
				.withLongOpt("angle")
				.withType(Number.class)
				.withDescription("angle increment, in degrees [default 5]")
				.create("a"));
		StringBuilder descr = new StringBuilder("Method to use. May be specified multiple times. [");
		RotationOrderMethod[] availableMethods = RotationOrderMethod.values();
		descr.append(availableMethods[0]);
		for(int i=1;i<availableMethods.length;i++) {
			descr.append(" | ");
			descr.append(availableMethods[i]);
		}
		descr.append(']');
		options.addOption(OptionBuilder.withArgName("method")
				.hasArg()
				.withLongOpt("method")
				.withType(RotationOrderMethod.class)
				.withDescription(descr.toString())
				.create('M') );
		options.addOption("o","output",true,"tab delimited output file. Angle " +
				"and Distance columns will be added. The string '%s' will be " +
				"expanded with the structure id.");
		options.addOption("d","display",false,"display jMol window with CeSymm alignment.");
		options.addOption(OptionBuilder.withLongOpt("nodisplay")
				.hasArg(false)
				.withDescription("don't display jMol window with CeSymm alignment.")
				.create() );
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


		String name;
		if(args.length == 0) {
			// default name
			name = "d1ijqa1";
			//		name = "1G6S";
			name = "1MER.A";
			//		name = "1MER";
			//		name = "1TIM.A";
			//		name = "d1h70a_";
			//name = "2YMS";

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
		} else if(args.length == 1) {
			name = args[0];
		} else {
			help.printHelp(usage, header, options, "");
			System.exit(1);
			return;
		}

		int maxorder = 8;
		if(cli.hasOption('x') ) {
			maxorder = Integer.parseInt(cli.getOptionValue('x'));
		}

		int angleIncr = 5; //degrees
		if(cli.hasOption('a') ) {
			angleIncr = Integer.parseInt(cli.getOptionValue('a'));
		}
		String outfile = null;
		if(cli.hasOption('o')) {
			outfile = String.format(cli.getOptionValue('o'),name);
		}

		List<OrderDetector> methods = new ArrayList<OrderDetector>();
		if(cli.hasOption('M')) {
			for(String method: cli.getOptionValues('M')) {
				RotationOrderMethod m = RotationOrderMethod.valueOf(method);
				methods.add(new RotationOrderDetector(maxorder,m) );
			}
		}
		//methods.add(new AngleOrderDetectorPlus(maxorder,100));
		//methods.add(new HybridOrderDetector(maxorder, Math.PI/16, false, .85));
		boolean displayAlignment = cli.hasOption('d') && ! cli.hasOption("nodisplay");

//		System.out.println("Name:" + name);
//		System.out.println("order:" + maxorder);
//		System.out.println("output:" + outfile);
//		for(RotationOrderDetector m: methods) {
//			System.out.println("method:"+m);
//		}
//		System.exit(0);

		// Done parsing arguments


		try {

			// Perform alignment to determine axis
			Atom[] ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
			Atom[] ca2 = StructureTools.cloneCAArray(ca1);
			CeSymm ce = new CeSymm();
			AFPChain alignment = ce.align(ca1, ca2);
			alignment.setName1(name);alignment.setName2(name);
			RotationAxis axis = new RotationAxis(alignment);

			// Output raw data
			if(outfile != null) {
				PrintStream out = null;
				try {
					if(outfile.equals("-")) {
						out = System.out;
					} else {
						out = new PrintStream(outfile);
					}
					out.println("Angle\tDistance");
					printSuperpositionDistance(ca1,axis,angleIncr*Calc.radiansPerDegree,out);
				} catch(FileNotFoundException e) {
					e.printStackTrace();
				} finally {
					if(out != null && out != System.out) {
						out.close();
					}
				}
			}
			// Display alignment
			if( displayAlignment ) {
				StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(alignment, ca1, ca2);
				jmol.evalString(axis.getJmolScript(ca1));
			}

			// Print orders
			for( OrderDetector detector: methods) {
				// Calculate order
				int order = detector.calculateOrder(alignment, ca1);
				System.out.format("%s\t%d%n",detector,order);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (StructureException e) {
			e.printStackTrace();
		} catch (OrderDetectionFailedException e) {
			e.printStackTrace();
		}

	}
}
