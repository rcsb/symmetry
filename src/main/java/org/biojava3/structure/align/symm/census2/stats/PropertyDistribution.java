package org.biojava3.structure.align.symm.census2.stats;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

public class PropertyDistribution {

	public static final String NEWLINE;
	private static final int MAX_FRACTION_DIGITS = 5;

	private static NumberFormat nf = new DecimalFormat();

	static final Logger logger = Logger.getLogger(BasicStats.class.getPackage().getName());

	static {
		NEWLINE = System.getProperty("line.separator");
	}
	static {
		BasicConfigurator.configure();
	}

	static {
		nf.setMaximumFractionDigits(MAX_FRACTION_DIGITS);
		nf.setMinimumFractionDigits(1);
	}

	public static void main(String[] args) {
		if (args.length != 2) {
			System.err.println("Usage: " + PropertyDistribution.class.getSimpleName() + " census-file significance-method");
			return;
		}
		File censusFile = new File(args[0]);
		Significance sig = SignificanceFactory.fromMethod(SignificanceFactory.class.getName(), args[1]);
		List<Property> properties = new ArrayList<Property>();
		properties.add(Property.identity());
		properties.add(Property.similarity());
		printBasicStats(System.out, censusFile, properties, sig);
	}

	public static void printBasicStats(PrintStream ps, File censusFile, List<Property> properties, Significance sig) {
		printBasicStats(new PrintWriter(ps, true), censusFile, properties, sig);
	}

	public static void printBasicStats(PrintWriter pw, File censusFile, List<Property> properties, Significance sig) {
		Grouping normalizingGrouping = Grouping.superfamily();
		Grouping reportGrouping = Grouping.fold();
		Results census;
		logger.info("Reading census file " + censusFile);
		try {
			census = Results.fromXML(censusFile);
		} catch (IOException e) {
			throw new RuntimeException("Could not load census file", e);
		}
		Results filtered = new Results();
		for (Result r : census.getData()) {
			if (sig.isSignificant(r)) filtered.add(r);
		}
		logger.info("Done reading census file");
		for (Property property : properties) {
			PropertyDistribution dist = new PropertyDistribution(normalizingGrouping, reportGrouping, filtered, property);
			pw.println(dist.toString());
			pw.println();
		}
	}

	private Grouping normalizingGrouping;
	private Grouping reportGrouping;
	private NavigableMap<String,List<Double>> map;
	private Property property;

	public PropertyDistribution(Grouping normalizingGrouping, Grouping reportGrouping, Results census, Property property) {

		this.property = property;
		this.normalizingGrouping = normalizingGrouping;
		this.reportGrouping = reportGrouping;
		map = new TreeMap<String,List<Double>>();
		Map<String,String> normalizingToReporting = new TreeMap<String,String>();

		// find values to normalize
		logger.info("Getting values of " + property.toString() + " from census.");
		Map<String,Double> values = new TreeMap<String,Double>();
		Map<String,Integer> counts = new TreeMap<String,Integer>();
		for (Result result : census.getData()) {
			String normalizingKey = normalizingGrouping.group(result);
			String reportKey = reportGrouping.group(result);
			try {
				double value = property.getProperty(result);
				StatUtils.plusD(values, normalizingKey, value);
				StatUtils.plus(counts, normalizingKey);
			} catch (PropertyUndefinedException e) {
				continue; // okay
			}
			normalizingToReporting.put(normalizingKey, reportKey);
		}

		// normalize (take mean)
		logger.info("Normalizing values over " + normalizingGrouping.toString());
		Map<String,Double> means = new TreeMap<String,Double>();
		for (Map.Entry<String,Double> entry : values.entrySet()) {
			final String normalizingKey = entry.getKey();
			final Double value = entry.getValue();
			final double fraction = value / (double) counts.get(normalizingKey);
			StatUtils.plusD(means, normalizingKey, fraction);
		}

		// okay, now record normalized stats
		logger.info("Reporting normalized values by " + reportGrouping.toString());
		for (Map.Entry<String,Double> entry : means.entrySet()) {
			final String normalizingKey = entry.getKey();
			final Double value = entry.getValue();
			final String reportKey = normalizingToReporting.get(normalizingKey);
			if (!map.containsKey(reportKey)) {
				map.put(reportKey, new ArrayList<Double>());
			}
			map.get(reportKey).add(value);
		}

	}

	public NavigableMap<String, List<Double>> getMap() {
		return map;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(property.toString() + ":" + NEWLINE);
		// for every fold, print a list of superfamily means
		for (Map.Entry<String,List<Double>> entry : map.entrySet()) {
			final String key = entry.getKey();
			final List<Double> values = entry.getValue();
			sb.append(key);
			for (double value : values) {
				sb.append("\t" + nf.format(value));
			}
			sb.append(NEWLINE);
		}
		return sb.toString();
	}

}
