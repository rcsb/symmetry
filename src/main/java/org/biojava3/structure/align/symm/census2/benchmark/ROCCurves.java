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
package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Generates and plots ROC curves for {@link Criterion Criteria} on benchmark {@link Case Cases}.
 * @author dmyerstu
 */
public class ROCCurves {

	static final Logger logger = Logger.getLogger(Case.class.getPackage().getName());

	static {
		BasicConfigurator.configure();
		logger.setLevel(Level.DEBUG);
	}

	private static final int DEFAULT_WIDTH = 1000;

	private static final int DEFAULT_HEIGHT = 800;

	public static void main(String[] args) {
		run(new File(args[0]), new File(args[1]));
	}

	public static void run(File input, File output) {
		ROCCurves rocs;
		try {
			List<Criterion<?>> criteria = new ArrayList<Criterion<?>>();
			criteria.add(Criterion.zScore(3.5f));
			criteria.add(Criterion.zScore(3.7f));
			criteria.add(Criterion.zScore(3.9f));
			criteria.add(Criterion.tmScore(0.3f));
			criteria.add(Criterion.tmScore(0.4f));
			criteria.add(Criterion.tmScore(0.5f));
			criteria.add(Criterion.inverseF(Criterion.screw(1.5f)));
			criteria.add(Criterion.random(0.5f));
//			criteria.add(Criterion.tmScore(0.6f));
//			Criterion<Float> c1 = Criterion.combineFF(Criterion.tmScore(0.3f), Criterion.zScore(3.5f), 2, 1);
//			criteria.add(c1);
//			criteria.add(Criterion.identity(0.15f));
//			criteria.add(Criterion.identity(0.1f));
//			criteria.add(Criterion.similarity(0.3f));
//			criteria.add(Criterion.epsilon((float) (Math.PI/16.0)));
//			criteria.add(Criterion.alignLength(100));
			rocs = new ROCCurves(input, criteria);
			rocs.graph(output);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private List<Criterion<?>> criteria;
	private Sample sample;

	public ROCCurves(File sampleFile, List<Criterion<?>> criteria) throws IOException {
		this(Sample.fromXML(sampleFile), criteria);
	}
	public ROCCurves(Sample sample, List<Criterion<?>> criteria) {
		List<Case> cases = new ArrayList<Case>(sample.size());
		for (Case c : sample.getData()) {
			if (c.getAlignment() == null || c.getAxis() == null) {
				logger.error("No alignment info for " + c.getScopId());
				continue;
			}
			cases.add(c);
			if (c.hasKnownSymmetry()) {
				nSymm++;
			} else {
				nAsymm++;
			}
		}
		sample.setData(cases);
		this.sample = sample;
		this.criteria = criteria;
	}

	private int nSymm;
	
	private int nAsymm;
	
	public SortedSet<Case> resort(final Criterion<?> criterion) {

		Comparator<Case> comparator = new Comparator<Case>() {
			Random random = new Random();
			@Override
			public int compare(Case o1, Case o2) {

				// this is the only case where we'll return 0
				if (o1.equals(o2)) return 0;

				// always put cases of known symmetry first
				try {
					if (criterion.hasSymmetry(o1.getResult()) && !criterion.hasSymmetry(o2.getResult())) return -1;
					if (!criterion.hasSymmetry(o1.getResult()) && criterion.hasSymmetry(o2.getResult())) return 1;
				} catch (NoncomputableCriterionException e) {
					throw new IllegalArgumentException(e);
				}

				// we don't want cases sorted
				return random.nextBoolean()? 1 : -1;

			}
		};

		SortedSet<Case> cases = new TreeSet<Case>(comparator);
		for (Case c : sample.getData()) {
			try {
				criterion.get(c.getResult());
			} catch (NoncomputableCriterionException e) {
				logger.warn("Can't compute " + criterion.getName() + " on " + c.getScopId());
				continue;
			}
			cases.add(c);
		}
		return cases;
	}

	public void printText(Criterion<?> criterion, String file) throws IOException {
		printText(criterion, new File(file));
	}
	public void printText(Criterion<?> criterion, File file) throws IOException {
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
		for (Case c : sample.getData()) {
			if (c.getResult() == null || c.getAlignment() == null) {
				pw.close();
				throw new DataIncompleteException(c.getScopId());
			}
			Number value;
			try {
				value = criterion.get(c.getResult());
			} catch (NoncomputableCriterionException e) {
				e.printStackTrace();
				continue;
			}
			pw.println(value + "\t" + (c.hasKnownSymmetry()? 1 : 0));
		}
		pw.close();
	}

	public void graph(String file) throws IOException {
		graph(new File(file), DEFAULT_WIDTH, DEFAULT_HEIGHT);
	}
	public void graph(String file, int width, int height) throws IOException {
		graph(new File(file), width, height);
	}
	public void graph(File file) throws IOException {
		graph(file, DEFAULT_WIDTH, DEFAULT_HEIGHT);
	}
	public void graph(File file, int width, int height) throws IOException {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (Criterion<?> criterion : criteria) {
			XYSeries series = new XYSeries(criterion.getName());
			int tp = 0, fp = 0;
			SortedSet<Case> cases = resort(criterion);
			logger.info("Adding series " + criterion.getName() + " with " + nSymm + " symmetric and " + nAsymm + " asymmetric:");
			for (Case c : cases) {
				double x, y;
				if (c.hasKnownSymmetry()) {
					tp++;
				} else {
					fp++;
				}
				y = (double) tp / (double) nSymm;
				x = (double) fp / (double) nAsymm;
				NumberFormat nf = new DecimalFormat();
				nf.setMaximumFractionDigits(3);
				if (c.hasKnownSymmetry()) {
					logger.debug(c.getScopId() + " (" + c.getKnownInfo() + ")" + " is symmetric: (" + nf.format(x) + ", " + nf.format(y) + ")");
				} else {
					logger.debug(c.getScopId() + " (" + c.getKnownInfo() + ")" + " is asymmetric: (" + nf.format(x) + ", " + nf.format(y) + ")");
				}
				series.add(x, y);
			}
			series.setDescription(criterion.getName());
			dataset.addSeries(series);
		}
		JFreeChart chart = ChartFactory.createXYLineChart("ROC", "FP", "TP", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYSplineRenderer renderer = new XYSplineRenderer();
		renderer.setShapesVisible(false);
		chart.getXYPlot().setRenderer(renderer);
		ChartUtilities.saveChartAsPNG(file, chart, width, height);
	}

}
