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
 * Created on 2013-03-10
 *
 */
package org.biojava.nbio.structure.align.symm.census3.stats;

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

import org.biojava.nbio.structure.align.symm.census3.CensusResult;
import org.biojava.nbio.structure.align.symm.census3.CensusResultList;
import org.biojava.nbio.structure.align.symm.census3.CensusSignificance;
import org.biojava.nbio.structure.align.symm.census3.CensusSignificanceFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Prints out the values of a certain {@link CensusResultProperty} in a {@link CensusResult}.
 * @author dmyersturnbull
 */
public class CensusPropertyDistribution {

	private final static Logger logger = LoggerFactory.getLogger(CensusPropertyDistribution.class);

	public static final String NEWLINE;
	private static final int MAX_FRACTION_DIGITS = 5;

	private static NumberFormat nf = new DecimalFormat();

	static {
		NEWLINE = System.getProperty("line.separator");
	}
	
	static {
		nf.setMaximumFractionDigits(MAX_FRACTION_DIGITS);
		nf.setMinimumFractionDigits(1);
	}

	public static void main(String[] args) {
		if (args.length != 2) {
			System.err.println("Usage: " + CensusPropertyDistribution.class.getSimpleName() + " census-file significance-method");
			return;
		}
		File censusFile = new File(args[0]);
		CensusSignificance sig = CensusSignificanceFactory.fromMethod(CensusSignificanceFactory.class.getName(), args[1]);
		List<CensusResultProperty<?>> properties = new ArrayList<CensusResultProperty<?>>();
		properties.add(CensusResultPropertyFactory.identity());
		printBasicStats(System.out, censusFile, properties, sig);
	}

	public static void printBasicStats(PrintStream ps, File censusFile, List<CensusResultProperty<?>> properties, CensusSignificance sig) {
		printBasicStats(new PrintWriter(ps, true), censusFile, properties, sig);
	}

	public static void printBasicStats(PrintWriter pw, File censusFile, List<CensusResultProperty<?>> properties, CensusSignificance sig) {
		StructureClassificationGrouping normalizingGrouping = StructureClassificationGrouping.superfamily();
		StructureClassificationGrouping reportGrouping = StructureClassificationGrouping.fold();
		CensusResultList census;
		logger.info("Reading census file " + censusFile);
		try {
			census = CensusResultList.fromXML(censusFile);
		} catch (IOException e) {
			throw new RuntimeException("Could not load census file", e);
		}
		CensusResultList filtered = new CensusResultList();
		for (CensusResult r : census.getEntries()) {
			if (sig.isSignificant(r)) filtered.add(r);
		}
		logger.info("Done reading census file");
		for (CensusResultProperty<?> property : properties) {
			CensusPropertyDistribution dist = new CensusPropertyDistribution(normalizingGrouping, reportGrouping, filtered, property);
			pw.println(dist.toString());
			pw.println();
		}
	}

	private StructureClassificationGrouping normalizingGrouping;
	private StructureClassificationGrouping reportGrouping;
	private NavigableMap<String,List<Double>> map;
	private CensusResultProperty<?> property;

	public CensusPropertyDistribution(StructureClassificationGrouping normalizingGrouping, StructureClassificationGrouping reportGrouping, CensusResultList census, CensusResultProperty property) {

		this.property = property;
		this.normalizingGrouping = normalizingGrouping;
		this.reportGrouping = reportGrouping;
		map = new TreeMap<String,List<Double>>();
		Map<String,String> normalizingToReporting = new TreeMap<String,String>();

		// find values to normalize
		logger.info("Getting values of " + property.toString() + " from census.");
		Map<String,Double> values = new TreeMap<String,Double>();
		Map<String,Integer> counts = new TreeMap<String,Integer>();
		for (CensusResult result : census.getEntries()) {
			String normalizingKey = normalizingGrouping.group(result);
			String reportKey = reportGrouping.group(result);
			try {
				double value = property.getProperty(result).doubleValue();
				CensusStatUtils.plusD(values, normalizingKey, value);
				CensusStatUtils.plus(counts, normalizingKey);
			} catch (PropertyUndefinedException e) {
				continue; // okay
			}
			normalizingToReporting.put(normalizingKey, reportKey);
		}

		// normalize (take mean)
		logger.info("Normalizing values over " + normalizingGrouping.toString());
		Map<String,Double> means = new TreeMap<String,Double>();
		double total = 0;
		for (Map.Entry<String,Double> entry : values.entrySet()) {
			final String normalizingKey = entry.getKey();
			final Double value = entry.getValue();
			final double fraction = value / (double) counts.get(normalizingKey);
			total += fraction;
			CensusStatUtils.plusD(means, normalizingKey, fraction);
		}

		System.err.println(total / means.size());
		
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
