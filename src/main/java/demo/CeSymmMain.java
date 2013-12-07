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
		final String usage = "[OPTIONS] [structures...]";
		final String header = "Determine the order for each structure, which may " +
				"be PDB IDs, SCOP domains, or file paths. If none are given, the " +
				"user will be prompted at startup.";

		Options options = new Options();
		options.addOption("h","help",false,"print help");
		options.addOption("j","jmol",false,"enable jMol display. [default]");
		options.addOption("J","nojmol",false,"disable jMol display.");
		options.addOption("a","alignment",false,"print alignment.");
		options.addOption("A","noalignment",false,"don't print alignment [default]");
		
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
		// Done parsing arguments


		for(String name: names) {
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


				// Order
				int symmNr = CeSymm.getSymmetryOrder(alignment);
				if(alignment.getTMScore()<.4) {
					symmNr = 1;
				}

				// Print result
				System.out.format("%s\tTMscore %f\tOrder %d%n",name,alignment.getTMScore(),symmNr);
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
}
