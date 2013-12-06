package demo;

import java.io.IOException;

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
		final String usage = "[OPTIONS] [structure]";
		final String header = "Determine the order for <structure>, which may " +
				"be a PDB ID, SCOP domain, or file path. If none is given, the " +
				"user will be prompted at startup.";

		Options options = new Options();
		options.addOption("h","help",false,"print help");
		options.addOption("J","nojmol",false,"disable jMol display.");
		options.addOption("j","jmol",false,"enable jMol display. [default]");
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
		} else if(args.length == 1) {
			name = args[0];
		} else {
			help.printHelp(usage, header, options, "");
			System.exit(1);
			return;
		}


		boolean displayAlignment = cli.hasOption('j') || !cli.hasOption('J');

		// Done parsing arguments


		try {

			// Perform alignment to determine axis
			Atom[] ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
			Atom[] ca2 = StructureTools.cloneCAArray(ca1);
			CeSymm ce = new CeSymm();
			AFPChain alignment = ce.align(ca1, ca2);
			RotationAxis axis = new RotationAxis(alignment);

			// Display alignment
			if( displayAlignment ) {
				StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(alignment, ca1, ca2);
				jmol.evalString(axis.getJmolScript(ca1));
			}

			// Print alignment
			System.out.println(alignment.toFatcat(ca1,ca2));

			// Print Order
			int symmNr = CeSymm.getSymmetryOrder(alignment);
			if(alignment.getTMScore()<.4) {
				symmNr = 1;
			}
			System.out.println("Symmetry order of: " + symmNr);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (StructureException e) {
			e.printStackTrace();
		}
	}
}
