package org.biojava3.structure.align.symm.census2.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

/**
 * Just makes the census fit for distribution by removing unnecessary information such as align-score and classification.
 * @author dmyerstu
 */
public class SimplifyCensusXML {

	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: " + SimplifyCensusXML.class.getSimpleName() + " census-input-file.xml census-output-file");
			return;
		}
		Results census = Results.fromXML(args[0]);
		File output = new File(args[1]);
		simplify(census);
		write(census, output);
	}

	public static void write(Results census, File output) throws IOException {
		BufferedWriter bw = null;
		File tmp = File.createTempFile(output.getName(), "tmp");
		try {
			bw = new BufferedWriter(new FileWriter(tmp));
			bw.write(census.toXML());
		} finally {
			if (bw != null) bw.close();
		}
		spacesToTabs(tmp, output, 4);
	}

	public static void simplify(Results census) {

		Significance sig = SignificanceFactory.forCeSymmOrd();

		for (Result result : census.getData()) {
			result.setClassification(null);
			result.setDescription(null);
			result.setFractionHelical(null);
			result.setIsSignificant(null);
			result.setSunId(null);
			result.setRank(null);
			result.getAlignment().setAlignScore(null);
			result.getAlignment().setCoverage(null);
			result.getAlignment().setzScore(null);
			result.getAlignment().setRmsd(null);
			result.getAlignment().setBlock2Length(null);
			result.getAlignment().setBlock1Length(null);

			// also change order
//			if (result.getOrder() == null || result.getOrder() < 2) {
//				if (result.getAxis() != null) {
//					int order = result.getAxis().guessOrder();
//					result.setOrder(order); // has a side effect of changing any -1s to 1
//				}
//			}

			result.setIsSignificant(sig.isSignificant(result));

		}

	}

	private static void spacesToTabs(File input, File output, int nSpaces) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(input));
		PrintWriter pw = new PrintWriter(output);
		String line = "";
		while ((line = br.readLine()) != null) {
			String trimmed = line.trim();
			int indent = (int) ((float) (line.length() - trimmed.length()) / (float) nSpaces);
			pw.println(repeat("\t", indent) + trimmed);
		}
		br.close();
		pw.close();
	}

	/**
	 * @return {@code s} repeated {@code n} times
	 */
	private static String repeat(String s, int n) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < n; i++) {
			sb.append(s);
		}
		return sb.toString();
	}

}
