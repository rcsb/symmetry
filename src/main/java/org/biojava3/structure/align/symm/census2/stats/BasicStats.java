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
package org.biojava3.structure.align.symm.census2.stats;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

/**
 * 
 * @author dmyerstu
 * 
 */
public class BasicStats {

	public static final String NEWLINE;
	private static final int MAX_FRACTION_DIGITS = 2;

	private static NumberFormat nf = new DecimalFormat();

	static final Logger logger = Logger.getLogger(BasicStats.class.getPackage().getName());

	private final Grouping grouping;

	private int nTotal;
	private int nSymm;
	private int nUnknown;
	private int nOrderUnknown;
	private int nRotational;

	private final NavigableMap<String, Integer> nTotalGroup = new TreeMap<String, Integer>();
	private final NavigableMap<String, Integer> nSymmGroup = new TreeMap<String, Integer>();
	private final NavigableMap<String, Integer> nUnknownGroup = new TreeMap<String, Integer>();
	private final NavigableMap<String, Integer> nRotGroup = new TreeMap<String, Integer>();
	private final NavigableMap<String, Integer> nRotUnknownGroup = new TreeMap<String, Integer>();

	private final NavigableMap<String, Integer> nTotalFold = new TreeMap<String, Integer>();
	private final NavigableMap<String, Integer> nSymmFold = new TreeMap<String, Integer>();
	private final NavigableMap<String, Integer> nRotFold = new TreeMap<String, Integer>();

	private final NavigableMap<String, Integer> nTotalClass = new TreeMap<String, Integer>();
	private final NavigableMap<String, Integer> nSymmClass = new TreeMap<String, Integer>();
	private final NavigableMap<String, Integer> nRotClass = new TreeMap<String, Integer>();

	private String text;

	protected Significance rotationallySymmetric = SignificanceFactory.getRotationallySymmetric();
	protected Significance sig = SignificanceFactory.getGenerallySymmetric();
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
		Grouping grouping = Grouping.byName(args[0]);
		File[] files = new File[args.length-1];
		for (int i = 1; i < args.length; i++) {
			files[i-1] = new File(args[i]);
		}
		printBasicStats(System.out, grouping, files);
	}

	public static void printBasicStats(PrintStream ps, Grouping grouping, File... censusFiles) {
		for (File file : censusFiles) {
			System.out.println(file.getName() + ":");
			printBasicStats(System.out, grouping, file);
			System.out.println();
		}
	}

	public static void printBasicStats(PrintStream ps, Grouping grouping, File censusFile) {
		BasicStats stats;
		try {
			stats = new BasicStats(grouping, censusFile);
		} catch (IOException e) {
			throw new RuntimeException("Could not load census file", e);
		}
		ps.println(stats.toString());
	}

	static void plus(Map<String, Double> map, String key, double value) {
		if (!map.containsKey(key)) map.put(key, 0.0);
		map.put(key, map.get(key) + value);
	}

	static void plus(Map<String, Integer> map, String key) {
		plus(map, key, 1);
	}

	static void plus(Map<String, Integer> map, String key, int value) {
		if (!map.containsKey(key)) map.put(key, 0);
		map.put(key, map.get(key) + value);
	}

	public BasicStats(Grouping grouping, File censusFile) throws IOException {
		this(grouping, Results.fromXML(censusFile));
	}

	public BasicStats(Grouping grouping, Results census) {

		/*
		 * First, collect stats per group:
		 * How many domains per group (e.g. superfamily) are symmetric?
		 */

		this.grouping = grouping;

		for (Result result : census.getData()) {

			final String group = grouping.group(result);

			boolean isSymm = false, isRot = false, isUnknown = false, isRotUnknown = false;

			nTotal++;

			// N symmetric and unknown symmetry
			try {
				if (sig.isSignificant(result) || sig == null && result.getIsSignificant()) {
					isSymm = true;
					nSymm++;
				}
			} catch (IllegalArgumentException e) {
				isUnknown = true;
				nUnknown++;
			}

			// N rotationally symmetric
			try {
				if (rotationallySymmetric.isSignificant(result)) {
					isRot = true;
					nRotational++;
				}
			} catch (IllegalArgumentException e) {
				isRotUnknown = true;
				if (isSymm) { // NOTE we assume can't find symmetry => can't find rotational symmetry
					nOrderUnknown++;
				}
			}

			// record stats per group
			plus(nTotalGroup, group);
			if (isSymm) plus(nSymmGroup, group);
			if (isRot) plus(nRotGroup, group);
			if (isUnknown) plus(nUnknownGroup, group);
			if (isRotUnknown) plus(nRotUnknownGroup, group);

			getAdditionalStats(result);

		}

		/*
		 * Now collect stats grouped by fold:
		 * How many groups (e.g. superfamilies) in the fold are symmetric?
		 */

		for (Map.Entry<String, Integer> entry : nTotalGroup.entrySet()) {

			final String group = entry.getKey();
			final String fold;
			try {
				String[] parts = group.split("\\.");
				fold = parts[0] + "." + parts[1]; // should always be fold
			} catch (ArrayIndexOutOfBoundsException e) {
				logger.error("Bad classification: " + group, e);
				continue;
			}

			// N total
			final int nTotalG = entry.getValue();
			plus(nTotalFold, fold);

			// % symmetric
			Integer nSymmG = nSymmGroup.get(group);
			if (nSymmG == null) nSymmG = 0;
			double fractionSymmG = (double) nSymmG / (double) nTotalG;
			if (fractionSymmG >= 0.5) {
				plus(nSymmFold, fold);
			}

			// % rotational
			Integer nRotG = nRotGroup.get(group);
			if (nRotG == null) nRotG = 0;
			double fractionRotG = (double) nRotG / (double) nTotalG;
			if (fractionRotG >= 0.5) {
				plus(nRotFold, fold);
			}

		}

		/*
		 * Finally, collect stats grouped by class:
		 * How many folds in each class are symmetric?
		 */

		for (Map.Entry<String, Integer> entry : nTotalFold.entrySet()) {

			final String fold = entry.getKey();
			final String clas;
			try {
				String[] parts = fold.split("\\.");
				clas = parts[0]; // should always be class
			} catch (ArrayIndexOutOfBoundsException e) {
				logger.error("Bad classification: " + fold, e);
				continue;
			}

			// N total
			final int nTotalG = entry.getValue();
			plus(nTotalClass, clas);

			// % symmetric
			Integer nSymmG = nSymmFold.get(fold);
			if (nSymmG == null) nSymmG = 0;
			double fractionSymmG = (double) nSymmG / (double) nTotalG;
			if (fractionSymmG >= 0.5) {
				plus(nSymmClass, clas);
			}

			// % rotational
			Integer nRotG = nRotFold.get(fold);
			if (nRotG == null) nRotG = 0;
			double fractionRotG = (double) nRotG / (double) nTotalG;
			if (fractionRotG >= 0.5) {
				plus(nRotClass, clas);
			}

		}

	}

	public final double getFractionOrderUnknown() {
		return (double) nOrderUnknown / (double) nTotal;
	}

	public final double getFractionRotational() {
		return (double) nRotational / (double) nTotal;
	}

	public final double getFractionSymmetric() {
		return (double) nSymm / (double) nTotal;
	}

	public final double getFractionUnknown() {
		return (double) nUnknown / (double) nTotal;
	}

	public final String getPercentageOrderUnknown() {
		return nf.format(getFractionOrderUnknown() * 100.0) + "%";
	}

	public final String getPercentageRotational() {
		return nf.format(getFractionRotational() * 100.0) + "%";
	}

	public final String getPercentageSymmetric() {
		return nf.format(getFractionSymmetric() * 100.0) + "%";
	}

	public final String getPercentageUnknown() {
		return nf.format(getFractionUnknown() * 100.0) + "%";
	}

	public final String getStatsText() {

		StringBuilder sb = new StringBuilder();

		// first print some basic per-domain stats
		// this is the only place where % unknown and % rotation unknown make sense
		sb.append("---PER DOMAIN---" + NEWLINE);
		sb.append("N total: " + nTotal + NEWLINE);
		sb.append("N symmetric: " + nSymm + NEWLINE);
		sb.append("% symmetric: " + getPercentageSymmetric() + NEWLINE);
		sb.append("N unknown: " + nUnknown + NEWLINE);
		sb.append("% unknown: " + getPercentageUnknown() + NEWLINE);
		sb.append("N rotational: " + nRotational + NEWLINE);
		sb.append("% rotational: " + getPercentageRotational() + NEWLINE);
		sb.append("N order unknown: " + nOrderUnknown + NEWLINE);
		sb.append("% order unknown: " + getPercentageOrderUnknown() + NEWLINE);
		sb.append(NEWLINE);

		/*
		 *  now we want to print grouped stats
		 */
		sb.append("---PER " + grouping.toString().toUpperCase() + "---" + NEWLINE);

		for (Map.Entry<String, Integer> entry : nTotalGroup.entrySet()) {

			final String group = entry.getKey();

			// N total
			final int nTotalG = entry.getValue();

			// % symmetric
			Integer nSymmG = nSymmGroup.get(group);
			if (nSymmG == null) nSymmG = 0;
			double fractionSymmG = (double) nSymmG / (double) nTotalG;

			// % rotational
			Integer nRotG = nRotGroup.get(group);
			if (nRotG == null) nRotG = 0;
			double fractionRotG = (double) nRotG / (double) nTotalG;

			// print for each group (e.g. superfamily)
			if (nSymmG > 0) {
				sb.append(group + ":" + NEWLINE);
				sb.append("\ttotal: " + nTotalG + NEWLINE);
				sb.append("\tsymmetric: " + fpc(fractionSymmG) + NEWLINE);
				sb.append("\trotational: " + fpc(fractionRotG) + NEWLINE);
			}

		}

		sb.append(NEWLINE);

		/* now report the % of symmetric superfamilies/groups per fold
		 * this is more important than the above, which is too detailed
		 * we also need to collect stats per class to print after
		 */
		sb.append("---PER " + grouping.toString() + " GROUPED BY FOLD---" + NEWLINE);

		for (Map.Entry<String, Integer> entry : nTotalFold.entrySet()) {

			final String fold = entry.getKey();

			// N total
			final int nTotalG = entry.getValue();

			// % symmetric
			Integer nSymmG = nSymmFold.get(fold);
			if (nSymmG == null) nSymmG = 0;
			double fractionSymmG = (double) nSymmG / (double) nTotalG;

			// % rotational
			Integer nRotG = nRotFold.get(fold);
			if (nRotG == null) nRotG = 0;
			double fractionRotG = (double) nRotG / (double) nTotalG;

			// print
			if (nSymmG > 0) {
				sb.append(fold + ":" + NEWLINE);
				sb.append("\ttotal: " + nTotalG + NEWLINE);
				sb.append("\tsymmetric: " + fpc(fractionSymmG) + NEWLINE);
				sb.append("\trotational: " + fpc(fractionRotG) + NEWLINE);
			}
		}

		sb.append(NEWLINE);

		/*
		 * now report the % of symmetric folds per class
		 */
		sb.append("---PER FOLD GROUPED BY CLASS---" + NEWLINE);

		for (Map.Entry<String, Integer> entry : nTotalClass.entrySet()) {

			final String clas = entry.getKey();

			// N total
			final int nTotalG = entry.getValue();

			// % symmetric
			Integer nSymmG = nSymmClass.get(clas);
			if (nSymmG == null) nSymmG = 0;
			double fractionSymmG = (double) nSymmG / (double) nTotalG;

			// % rotational
			Integer nRotG = nRotClass.get(clas);
			if (nRotG == null) nRotG = 0;
			double fractionRotG = (double) nRotG / (double) nTotalG;

			// print
			sb.append(clas + ":" + NEWLINE);
			sb.append("\ttotal: " + nTotalG + NEWLINE);
			sb.append("\tsymmetric: " + fpc(fractionSymmG) + NEWLINE);
			sb.append("\trotational: " + fpc(fractionRotG) + NEWLINE);

		}

		return sb.toString();

	}

	@Override
	public String toString() {
		if (text == null) text = getStatsText();
		return text;
	}

	private final String fpc(double p) {
		return nf.format(p * 100.0) + "%";
	}

	protected void getAdditionalStats(Result result) {
	}

}
