package org.biojava3.structure.align.symm.benchmark;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.biojava3.structure.align.symm.census3.CensusResultList;


/**
 * Gets a Results object (see symmetry project) back from a {@link Sample} object.
 * @author dmyerstu
 */
public class ResultsExtractor {

	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: " + ResultsExtractor.class.getSimpleName() + " input-benchmark-file output-results-file");
			return;
		}
		extractResults(new File(args[0]), new File(args[1]));
	}

	private static void extractResults(File input, File output) throws IOException {
		Sample sample = Sample.fromXML(input);
		CensusResultList results = extractResults(sample);
		String xml = results.toXML();
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		bw.write(xml);
		bw.close();
	}

	private static CensusResultList extractResults(Sample sample) {
		CensusResultList results = new CensusResultList();
		for (Case c : sample.getData()) {
			results.add(c.getResult());
		}
		return results;
	}

}
