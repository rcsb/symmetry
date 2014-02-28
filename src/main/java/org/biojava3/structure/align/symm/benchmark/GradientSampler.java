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
 * Created on 2013-03-28
 *
 */
package org.biojava3.structure.align.symm.benchmark;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.benchmark.comparison.AccuracyFinder;
import org.biojava3.structure.align.symm.benchmark.comparison.Criterion;
import org.biojava3.structure.align.symm.benchmark.comparison.ROCCurves;
import org.biojava3.structure.align.symm.census2.Census.AlgorithmGiver;
import org.biojava3.structure.align.symm.census2.NamesCensus;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;

/**
 * Runs CE-Symm with different gradient penalties on the main diagonal, and benchmark the results.
 * @author dmyerstu
 * @deprecated We found using a gradient to do very little
 */
public class GradientSampler {

	private static final Logger logger = LogManager.getLogger(GradientSampler.class.getName());

	private double[][] gradients;

	private double[] tp;
	private double[] fp;

	public GradientSampler(double[][] gradients) {
		this.gradients = gradients;
	}

	public GradientSampler(int order, double start, double stop, double step) {
		List<double[]> list = new ArrayList<double[]>();
		for (double x = start; x <= stop; x += step) {
			double[] gradient = new double[order];
			gradient[0] = x;
			list.add(gradient);
		}
		this.gradients = new double[list.size()][];
		for (int i = 0; i < list.size(); i++) {
			this.gradients[i] = list.get(i);
		}
	}

	public String gradToString(int x) {
		StringBuilder sb = new StringBuilder(x + "=(");
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(0);
		for (int i = 0; i < gradients[x].length; i++) {
			sb.append(nf.format(gradients[x][i]));
			if (i < gradients[x].length-1) sb.append(' ');
		}
		sb.append(")");
		return sb.toString();
	}
	
	public void sample(String pdbDir, String dir, File lineByLine, File ordersFile) {
		if (!dir.endsWith("/")) dir += "/";
		final int n = gradients.length;
		this.tp = new double[n];
		this.fp = new double[n];
		List<Criterion> criteria = new ArrayList<Criterion>();
		criteria.add(Criterion.combine(Criterion.tmScore(), Criterion.hasOrder(1), 1, 1));
		for (int i = 0; i < n; i++) {
			System.out.println("SAMPING: " + i);
			final int j = i;
			AlgorithmGiver algorithm = new AlgorithmGiver() {
				@Override
				public StructureAlignment getAlgorithm() {
					CeSymm ce = new CeSymm();
//					ce.setGradientPolyCoeff(gradients[j]);
					return ce;
				}
			};
			File censusFile = new File(dir + gradToString(i) + "_benchmark_stub.xml");
			NamesCensus.buildDefault(censusFile, lineByLine, algorithm);
			File benchmarkFile = new File(dir + i + "_benchmark.xml");
			try {
				SampleBuilder.buildSample(censusFile, benchmarkFile, ordersFile);
			} catch (IOException e) {
				throw new RuntimeException("Couldn't build benchmark", e);
			}
			try {
				File accFile = new File(dir + i + "_acc.txt");
				AccuracyFinder finder = new AccuracyFinder(benchmarkFile, CensusSignificanceFactory.forCeSymmOrd());
				PrintWriter accPw = new PrintWriter(new BufferedWriter(new FileWriter(accFile)));
				accPw.println(finder);
				accPw.close();
				tp[i] = finder.getTp();
				fp[i] = finder.getFp();
			} catch (IOException e) {
				throw new RuntimeException("Couldn't load benchmark", e);
			}
			try {
				File rocFile = new File(dir + i + "_roc.png");
				File rocTextFile = new File(dir + i + "_roc.txt");
				ROCCurves roc = new ROCCurves(benchmarkFile, criteria);
				roc.graph(rocFile);
				PrintWriter rocPw = new PrintWriter(new BufferedWriter(new FileWriter(rocTextFile)));
				roc.printText(rocPw);
				rocPw.close();
			} catch (IOException e) {
				logger.error("Couldn't make ROC curves for " + gradToString(i), e);
			}
		}
	}

	public double[] getMaxTp() {
		double max = Double.MIN_VALUE;
		int n = -1;
		for (int i = 0; i < gradients.length; i++) {
			if (tp[i] > max) {
				max = tp[i];
				n = i;
			}
		}
		return gradients[n];
	}

	public double[] getMinFp() {
		double min = Double.MAX_VALUE;
		int n = -1;
		for (int i = 0; i < gradients.length; i++) {
			if (fp[i] < min) {
				min = fp[i];
				n = i;
			}
		}
		return gradients[n];
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("gradient\ttp\tfp\n");
		for (int i = 0; i < gradients.length; i++) {
			sb.append(gradients[i] + "\t" + tp[i] + "\t" + fp[i] + "\n");
		}
		return sb.toString();
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String pdbDir = args[0];
		String mainDir = args[1];
		File lineByLine = new File(args[2]);
		File ordersFile = new File(args[3]);
		if (!mainDir.endsWith("/")) mainDir += "/";
		new File(mainDir).mkdirs();
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(mainDir + "results")));
		for (int i = 1; i < 3; i++) {
			try {
				GradientSampler sampler = new GradientSampler(i, 10, 1000, 100);
				String subDir = mainDir + i + "/";
				new File(subDir).mkdirs();
				sampler.sample(pdbDir, subDir, lineByLine, ordersFile);
				pw.println("FOR " + i);
				pw.println(sampler);
				pw.println("MAX(TP): " + sampler.getMaxTp());
				pw.println("MIN(FP): " + sampler.getMinFp());
				pw.println();
			} catch (RuntimeException e) {
				pw.println(e);
				e.printStackTrace(pw);
				logger.error("Couldn't run benchmark " + i, e);
			}
		}
		pw.close();
	}

}
