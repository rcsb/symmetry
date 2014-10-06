package org.biojava3.structure.align.symm.census3.stats;

import java.io.File;
import java.io.IOException;
import java.util.Formatter;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusResultList;
import org.biojava3.structure.align.symm.census3.CensusSignificance;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;

/**
 * Extremely simple statistics using the default {@link Significance} objects for CE-Symm results. Simply takes a list
 * of files and prints the percentage symmetry in each file, without any normalization (that is, just per domain).
 * 
 * @author dmyersturnbull
 */
public class NonNormalizedCensusStats {

	private static final Logger logger = LogManager.getLogger(NonNormalizedCensusStats.class.getName());

	private static CensusSignificance ord = CensusSignificanceFactory.forCeSymmOrd();

	private static CensusSignificance tm = CensusSignificanceFactory.forCeSymmTm();

	private Map<String, Integer> counts = new HashMap<String, Integer>();
	private Map<String, Double> ordStats = new HashMap<String, Double>();
	private Map<String, Double> tmStats = new HashMap<String, Double>();

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			System.err.println("Usage: " + NonNormalizedCensusStats.class.getSimpleName() + " file-1 [file-2 ...]");
		}
		File[] files = new File[args.length];
		for (int i = 0; i < args.length; i++)
			files[i] = new File(args[i]);
		NonNormalizedCensusStats stats = new NonNormalizedCensusStats(files);
		System.out.println(stats);
	}

	public NonNormalizedCensusStats(File... files) throws IOException {
		for (File file : files) {
			CensusResultList results = CensusResultList.fromXML(file);
			int nTm = 0, nOrd = 0;
			for (CensusResult result : results.getEntries()) {
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
		sb.append("file-name\tord\tTM\tN" + CensusStatUtils.NEWLINE);
		for (String file : tmStats.keySet()) {
			// formatter.format("%-8s %8.2f%% %8.2f%% %n", file, tmStats.get(file).doubleValue()*100.0,
			// ordStats.get(file).doubleValue()*100.0, counts.get(file).intValue());
			sb.append(file + "\t" + CensusStatUtils.formatP(ordStats.get(file)) + "\t" + CensusStatUtils.formatP(tmStats.get(file))
					+ "\t" + counts.get(file) + CensusStatUtils.NEWLINE);
		}
		formatter.close();
		return sb.toString();
	}

}
