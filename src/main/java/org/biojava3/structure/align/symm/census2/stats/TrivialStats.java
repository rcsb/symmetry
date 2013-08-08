package org.biojava3.structure.align.symm.census2.stats;

import java.io.File;
import java.io.IOException;
import java.util.Formatter;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

/**
 * Extremely simple statistics using the default {@link Significance} objects for CE-Symm results. Simply takes a list
 * of files and prints the percentage symmetry in each file, without any normalization (that is, just per domain).
 * 
 * @author dmyerstu
 */
public class TrivialStats {

	private static final Logger logger = LogManager.getLogger(TrivialStats.class.getPackage().getName());

	private static Significance ord = SignificanceFactory.rotationallySymmetricSmart();

	private static Significance tm = SignificanceFactory.generallySymmetric();

	private Map<String, Integer> counts = new HashMap<String, Integer>();
	private Map<String, Double> ordStats = new HashMap<String, Double>();
	private Map<String, Double> tmStats = new HashMap<String, Double>();

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			System.err.println("Usage: " + TrivialStats.class.getSimpleName() + " file-1 [file-2 ...]");
		}
		File[] files = new File[args.length];
		for (int i = 0; i < args.length; i++)
			files[i] = new File(args[i]);
		TrivialStats stats = new TrivialStats(files);
		System.out.println(stats);
	}

	public TrivialStats(File... files) throws IOException {
		for (File file : files) {
			Results results = Results.fromXML(file);
			int nTm = 0, nOrd = 0;
			for (Result result : results.getData()) {
				if (tm.isSignificant(result)) nTm++;
				if (ord.isSignificant(result)) nOrd++;
			}
			tmStats.put(file.getName(), (double) nTm / (double) results.size());
			ordStats.put(file.getName(), (double) nOrd / (double) results.size());
			counts.put(file.getName(), results.size());
		}
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		Formatter formatter = new Formatter(sb);
		sb.append("file-name\tord\tTM\tN" + StatUtils.NEWLINE);
		for (String file : tmStats.keySet()) {
			// formatter.format("%-8s %8.2f%% %8.2f%% %n", file, tmStats.get(file).doubleValue()*100.0,
			// ordStats.get(file).doubleValue()*100.0, counts.get(file).intValue());
			sb.append(file + "\t" + StatUtils.formatP(ordStats.get(file)) + "\t" + StatUtils.formatP(tmStats.get(file))
					+ "\t" + counts.get(file) + StatUtils.NEWLINE);
		}
		formatter.close();
		return sb.toString();
	}

}
