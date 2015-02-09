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
 * Created on 2013-02-22
 *
 */
package org.biojava.nbio.structure.align.symm.benchmark.comparison;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.nbio.structure.align.symm.benchmark.Case;
import org.biojava.nbio.structure.align.symm.benchmark.Sample;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Generates and plots ROC curves for {@link Criterion Criteria} on a benchmark {@link Sample}.
 * 
 * @author dmyerstu
 */
public class ROCCurves {

	private static final Logger logger = LogManager.getLogger(ROCCurves.class.getName());

	private static final int DEFAULT_HEIGHT = 1600;

	private static final int DEFAULT_WIDTH = 1600;

	private List<Criterion> criteria;

	private XYSeriesCollection dataset;

	private Sample sample;

	/**
	 * 
	 * @param args
	 * <ol>
	 * <li>An input XML file (see {@link Sample}) of CE-Symm benchmark data</li>
	 * <li>The path to which to output the ROC curves for CE-Symm</li>
	 * <li>An input XML file (see {@link Sample}) of SymD benchmark data</li>
	 * <li>The path to which to output the ROC curves for SymD</li>
	 * </ol>
	 * 
	 */
	public static void main(String[] args) {
		if (args.length < 2 || args.length > 5) {
			System.err.println("Usage: " + ROCCurves.class.getSimpleName() + " input-cesymm-benchmark-file output-cesymm-file [input-symd-benchmark-file output-symd-file [1.5b]]");
			return;
		}
		runForCeSymm(new File(args[0]), new File(args[1]));
		boolean updated = false;
		if (args.length > 4) {
			if (args[4].equalsIgnoreCase("1.5b") || args[4].equalsIgnoreCase("true") || args[4].equalsIgnoreCase("updated")) updated = true;
		}
		if (args.length > 2) {
			runForSymD(new File(args[2]), new File(args[3]), updated);
		}
	}

	/**
	 * Does two things:
	 * <ol>
	 * <li>Prints an ROC curve of the CE-Symm data to {@code ceSymmOutput}</li>
	 * <li>Prints a list of data points for CE-Symm to standard output (for graphing in a spreadsheet)</li>
	 * </ol>
	 * @param input A benchmark XML file (see {@link Sample}) containing data from CE-Symm
	 * @param ceSymmOutput
	 * @param symdOutput
	 */
	public static void runForCeSymm(File input, File output) {
		try {

			// ROC curves
			List<Criterion> ceSymmCriteria = new ArrayList<Criterion>();
			ceSymmCriteria.add(Criterion.tmScore());
			ceSymmCriteria.add(Criterion.combine(Criterion.tmScore(), Criterion.hasOrderLiberal(1f), 1, 1));
			ROCCurves ceSymmRocs = new ROCCurves(input, ceSymmCriteria);
			ceSymmRocs.graph(output);
			ceSymmRocs.printMatrices(System.out);
			// print text
			ceSymmRocs.printText();

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Does two things:
	 * <ol>
	 * <li>Prints an ROC curve of the SymD data to {@code ceSymmOutput}</li>
	 * <li>Prints a list of data points for SymD to standard output (for graphing in a spreadsheet)</li>
	 * </ol>
	 * @param input A benchmark XML file (see {@link Sample}) containing data from SymD
	 * @param ceSymmOutput
	 * @param symdOutput
	 */
	public static void runForSymD(File input, File output, boolean isUpdated) {
//		try {
//
//			// ROC curves
//			List<Criterion> symdCriteria = new ArrayList<Criterion>();
//			symdCriteria.add(Criterion.symdZScore());
//			symdCriteria.add(Criterion.tScore());
//			if (isUpdated) symdCriteria.add(Criterion.symDTMScore());
//			symdCriteria.add(Criterion.tmScore());
//			ROCCurves symdRocs = new ROCCurves(input, symdCriteria);
//			symdRocs.graph(output);
//
//			// print text
//			symdRocs.printMatrices(System.out);
//			symdRocs.printText();
//
//		} catch (IOException e) {
//			throw new RuntimeException(e);
//		}
	}

	/**
	 * @see #ROCCurves(Sample, List)
	 * @param sampleFile
	 * @param criteria
	 * @throws IOException
	 */
	public ROCCurves(File sampleFile, List<Criterion> criteria) throws IOException {
		this(Sample.fromXML(sampleFile), criteria);
	}

	/**
	 * Initializes a new list of ROC curve data with the specified criteria. If a {@link Case} in the {@link Sample} has an error, logs the error but does not stop.
	 * @param sample
	 * @param criteria
	 */
	public ROCCurves(Sample sample, List<Criterion> criteria) {
		List<Case> cases = new ArrayList<Case>(sample.size());
		for (Case c : sample.getData()) {
			try {
				c.getKnownInfo().hasRotationalSymmetry();
				cases.add(c); // only add if the known/benchmark result is known
			} catch (RuntimeException e) {
				logger.fatal("Encountered an error on " + c.getScopId(), e);
				throw e; // just be safe
			}
		}
		logger.info("Found " + sample.size() + " cases");
		sample.setData(cases);
		this.sample = sample;
		this.criteria = criteria;
	}

	public void graph(File file) throws IOException {
		graph(file, DEFAULT_WIDTH, DEFAULT_HEIGHT);
	}

	public void graph(File file, int width, int height) throws IOException {
		if (dataset == null) dataset = getRocPoints(); // lazy initialization
		JFreeChart chart = ChartFactory.createXYLineChart("ROC", "FP", "TP", dataset, PlotOrientation.VERTICAL, true,
				false, false);
		ChartUtilities.saveChartAsPNG(file, chart, width, height);
	}

	public void graph(String file) throws IOException {
		graph(new File(file), DEFAULT_WIDTH, DEFAULT_HEIGHT);
	}

	public void graph(String file, int width, int height) throws IOException {
		graph(new File(file), width, height);
	}

	public void printText() throws IOException {
		printText(System.out);
	}

	public void printText(PrintStream ps) {
		if (dataset == null) dataset = getRocPoints(); // lazy initialization

		for (int i = 0; i < criteria.size(); i++) {

			// print the name of the criterion
			final Criterion criterion = criteria.get(i);
			final XYSeries series = dataset.getSeries(i);
			ps.println(criterion);

			// print the X points
			for (int j = 0; j < series.getItemCount(); j++) {
				final double x = series.getX(j).doubleValue();
				ps.print(x + "\t");
			}
			ps.println();

			// print the Y points
			for (int j = 0; j < series.getItemCount(); j++) {
				final double y = series.getY(j).doubleValue();
				ps.print(y + "\t");
			}
			ps.println("\n");

		}
	}

	public void printMatrices(PrintStream ps) {
		List<Case> positive = new ArrayList<Case>();
		for (Case c : sample.getData()) {
			if (c.getKnownInfo().hasRotationalSymmetry()) {
				positive.add(c);
			}
		}
		List<Case> negative = new ArrayList<Case>();
		for (Case c : sample.getData()) {
			if (!c.getKnownInfo().hasRotationalSymmetry()) {
				negative.add(c);
			}
		}
		ps.println("------------------POSITIVE------------------");
		for (Case c : positive) {
			Iterator<Criterion> iter = criteria.iterator();
			while (iter.hasNext()) {
				double value = -1;
				try {
					value = iter.next().get(c.getResult());
				} catch (NoncomputableCriterionException e) {
					e.printStackTrace();
				}
				ps.print(value);
				if (iter.hasNext()) ps.print("\t");
			}
			ps.println();
		}
		ps.println(); ps.println();
		ps.println("------------------NEGATIVE------------------");
		for (Case c : negative) {
			Iterator<Criterion> iter = criteria.iterator();
			while (iter.hasNext()) {
				double value = -1;
				try {
					value = iter.next().get(c.getResult());
				} catch (NoncomputableCriterionException e) {
					e.printStackTrace();
				}
				ps.print(value);
				if (iter.hasNext()) ps.print("\t");
			}
			ps.println();
		}
	}

	/**
	 * For each {@link Criterion}, prints a line of X coordinates followed by a line of Y coordinates.
	 * @param ps
	 */
	public void printText(PrintWriter ps) {

		if (dataset == null) dataset = getRocPoints(); // lazy initialization

		for (int i = 0; i < criteria.size(); i++) {

			// print the name of the criterion
			final Criterion criterion = criteria.get(i);
			final XYSeries series = dataset.getSeries(i);
			ps.println(criterion);

			// print the X points
			for (int j = 0; j < series.getItemCount(); j++) {
				final double x = series.getX(j).doubleValue();
				ps.print(x + "\t");
			}
			ps.println();

			// print the Y points
			for (int j = 0; j < series.getItemCount(); j++) {
				final double y = series.getY(j).doubleValue();
				ps.print(y + "\t");
			}
			ps.println("\n");

		}
	}

	private List<Case> resort(final Criterion criterion) {

		// we want to sort by the scoring criterion
		Comparator<Case> comparator = new Comparator<Case>() {
			Random random = new Random();

			@Override
			public int compare(Case o1, Case o2) {

				// this is the only case where we'll return 0
				if (o1.equals(o2)) return 0;

				try {
					double c1 = criterion.get(o1.getResult());
					double c2 = criterion.get(o2.getResult());
					if (c1 < c2) return 1;
					if (c1 > c2) return -1;
				} catch (NoncomputableCriterionException e) {
					throw new IllegalArgumentException(e);
				}

				// select randomly when there's a tie
				return random.nextBoolean() ? 1 : -1;

			}
		};

		SortedSet<Case> cases = new TreeSet<Case>(comparator);
		for (Case c : sample.getData()) {
			try {
				criterion.get(c.getResult());
			} catch (NoncomputableCriterionException e) { // this shouldn't happen since we've already initialized
				logger.error("Can't compute " + criterion.getName() + " on " + c.getScopId());
				continue;
			}
			cases.add(c);
		}

		// use a list to decrease overall time-complexity
		List<Case> caseList = new ArrayList<Case>(cases.size());
		caseList.addAll(cases);
		return caseList;
	}

	public XYSeriesCollection getRocPoints() {

		if (dataset != null) return dataset;

		XYSeriesCollection dataset = new XYSeriesCollection();

		for (Criterion criterion : criteria) {

			// get the cases sorted by the current criterion
			// this should contain only computable values
			List<Case> cases = resort(criterion);
			XYSeries series = new XYSeries(criterion.getName());

			/*
			 * Instead of creating a new subsample each time,
			 * start with the minimal (size=0) subsample, then incrementally add cases to the subsample,
			 * and increment true positives (tp; y) or false positives (fp; x) each iteration
			 * This is equivalent to creating a totally new subsample each time,
			 * but it's O(nlogn) + O(n) rather than O(nlogn) + O(n^2).
			 */
			int tp = 0, fp = 0;

			/*
			 * We need updated statistics for the number of positive and negative known values.
			 * This is because data points may not exist for this criterion (or we may have new data points).
			 */
			int nSymm = 0, nAsymm = 0;
			for (Case c: cases) {
				if (c.getKnownInfo().hasRotationalSymmetry()) {
					nSymm++;
				} else {
					nAsymm++;
				}
			}

			for (Case c : cases) {

				// get the new x and y points
				double x, y;
				if (c.getKnownInfo().hasRotationalSymmetry()) {
					tp++;
				} else {
					fp++;
				}
				y = (double) tp / (double) nSymm;
				x = (double) fp / (double) nAsymm;

				// now just add the data point
				NumberFormat nf = new DecimalFormat();
				nf.setMaximumFractionDigits(3);
				if (c.getKnownInfo().hasRotationalSymmetry()) {
					logger.debug(c.getScopId() + " (" + c.getKnownInfo() + ")" + " is symmetric: (" + nf.format(x)
							+ ", " + nf.format(y) + ")");
				} else {
					logger.debug(c.getScopId() + " (" + c.getKnownInfo() + ")" + " is asymmetric: (" + nf.format(x)
							+ ", " + nf.format(y) + ")");
				}
				series.add(x, y);

			}

			logger.info("Adding series " + criterion.getName() + " with " + nSymm + " symmetric and " + nAsymm
					+ " asymmetric:");

			// add the series to the collection
			series.setDescription(criterion.getName());
			dataset.addSeries(series);

		}

		return dataset;
	}

}
