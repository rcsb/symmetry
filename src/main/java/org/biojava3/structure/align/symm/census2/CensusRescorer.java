/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2013-03-25
 *
 */
package org.biojava3.structure.align.symm.census2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * A class to update a census ({@link Results}) XML file with new score cutoffs.
 * @author dmyerstu
 */
public class CensusRescorer {

	private Significance sig;

	public CensusRescorer(Significance sig) {
		this.sig = sig;
	}

	public void rescore(Results census) {
		for (Result result : census.getData()) {
			result.setIsSignificant(sig.isSignificant(result));
		}
	}

	public static void convert(Significance sig, File input, File output) throws IOException {
		CensusRescorer rescorer = new CensusRescorer(sig);
		Results census = Results.fromXML(input);
		rescorer.rescore(census);
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		bw.write(census.toXML());
		bw.close();
	}

	/**
	 * @param args
	 * <ol>
	 * <li>The name of the method in SignificanceFactory to use for Significance</li>
	 * <li>The XML input file</li>
	 * <li>The path of the output file to write</li>
	 * </ol>
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length > 3 || args.length < 2) {
			System.out.println("Usage: CensusRescorer intput-file output-file [significance-method]");
			return;
		}
		File input = new File(args[0]);
		File output = new File(args[1]);
		Significance sig = SignificanceFactory.getGenerallySymmetric();
		if (args.length == 3) {
			String sigMethod = args[2];
			sig = SignificanceFactory.fromMethod(null, sigMethod);
		}
		convert(sig, input, output);
	}

}
