package org.biojava3.structure.align.symm.census2.stats;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;

public class SmartStats {

	static final Logger logger = Logger.getLogger(BasicStats.class.getPackage().getName());

	static {
		BasicConfigurator.configure();
	}

	private static double DEFAULT_TM_SCORE_THRESHOLD = 0.4;

	private static double DEFAULT_ORDER_THRESHOLD = 0.5;

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 3) {
			System.err.println("Usage: SmartStats census-file [tm-score-threshold] [order-threshold]");
			return;
		}
		File census = new File(args[0]);
		double tmScoreTreshold = DEFAULT_TM_SCORE_THRESHOLD;
		double orderTreshold = DEFAULT_ORDER_THRESHOLD;
		if (args.length > 1) tmScoreTreshold = Double.parseDouble(args[1]);
		if (args.length > 2) orderTreshold = Double.parseDouble(args[2]);
		SmartStats stats = new SmartStats(census, tmScoreTreshold, orderTreshold);
		System.out.println(stats);
	}
	
	private Property property;

	public Property getProperty() {
		return property;
	}
	public SmartStats(File census, double tmScoreTreshold, double orderTreshold) throws IOException {
		this(Property.tmScore(), Results.fromXML(census), tmScoreTreshold, orderTreshold);
	}
	public SmartStats(Property property, Results census, double tmScoreTreshold, double orderTreshold) {
		
		this.property = property;
		
		final Grouping sfGrouping = Grouping.superfamily();
		Map<String,Double> tmInSfs = new TreeMap<String,Double>();
		Map<String,Integer> nInSfs = new TreeMap<String,Integer>();
		Map<String,Integer> nWithOrderInSfs = new TreeMap<String,Integer>();
		
		for (Result result : census.getData()) {
			final String sf = sfGrouping.group(result);
			double tmScore = 0;
			int order = 0;
			StatUtils.plus(nInSfs, sf);
			final String fold;
			final String clas;
			try {
				String[] parts = sf.split("\\.");
				fold = parts[0] + "." + parts[1];
				clas = parts[0];
			} catch (ArrayIndexOutOfBoundsException e) {
				logger.error("Bad classification: " + sf, e);
				continue;
			}
			try {
				tmScore = property.getProperty(result);
			} catch (PropertyUndefinedException e) {
				// okay, just leave as 0
			}
			try {
				order = (int) Property.hasOrder().getProperty(result);
			} catch (PropertyUndefinedException e) {
				// okay, just leave as 0
			}
			StatUtils.plus(nDomainsInFolds, fold);
			StatUtils.plus(nDomainsInFolds, clas);
			StatUtils.plus(nDomainsInFolds, "overall");
			StatUtils.plusD(tmInSfs, sf, tmScore);
			StatUtils.plus(nWithOrderInSfs, sf, order);
		}
		
		
		for (Map.Entry<String, Integer> entry : nInSfs.entrySet()) {
			final String sf = entry.getKey();

			final String fold;
			final String clas;
			try {
				String[] parts = sf.split("\\.");
				fold = parts[0] + "." + parts[1];
				clas = parts[0];
			} catch (ArrayIndexOutOfBoundsException e) {
				logger.error("Bad classification: " + sf, e);
				continue;
			}
			final int nInSf = entry.getValue();
			final double fracWithOrder = nWithOrderInSfs.get(sf).doubleValue() / (double) nInSf;
			final double meanTmScore = tmInSfs.get(sf) / (double) nInSf;
			if (fracWithOrder >= orderTreshold && meanTmScore >= tmScoreTreshold) {
				StatUtils.plus(nSymmInFolds, fold);
				StatUtils.plus(nSymmInFolds, clas);
				StatUtils.plus(nSymmInFolds, "overall");
			}
			StatUtils.plus(nInFoldsTotal, fold);
			StatUtils.plus(nInFoldsTotal, clas);
			StatUtils.plus(nInFoldsTotal, "overall");
		}

	}

	Map<String,Integer> nDomainsInFolds = new TreeMap<String,Integer>();
	Map<String,Integer> nInFoldsTotal = new TreeMap<String,Integer>();
	Map<String,Integer> nSymmInFolds = new TreeMap<String,Integer>();
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String, Integer> entry : nInFoldsTotal.entrySet()) {
			final String fold = entry.getKey();
			final int total = entry.getValue();
			double fractionSymm;
			if (nSymmInFolds.get(fold) != null) {
				fractionSymm = (double) nSymmInFolds.get(fold) / (double) total;
			} else {
				fractionSymm = 0;
			}
			if (nDomainsInFolds.get(fold) >= 10 && fractionSymm >= 0.3 || fold.length() < 3 || fold.equals("overall")) {
				sb.append(fold + "\t" + nInFoldsTotal.get(fold) + "\t" + nDomainsInFolds.get(fold) + "\t" + StatUtils.formatP(fractionSymm) + StatUtils.NEWLINE);
			}
		}
		return sb.toString();
	}
	
}
