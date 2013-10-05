package org.biojava3.structure.align.symm.benchmark.comparison.order;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.benchmark.SampleBuilder;
import org.biojava3.structure.align.symm.census2.Census.AlgorithmGiver;
import org.biojava3.structure.align.symm.census2.NamesCensus;

/**
 * Finds optimal minMetricChange.
 * @author dmyersturnbull
 */
public class ThresholdTester {

	public static void main(String[] args) throws FileNotFoundException {
		run(0.1f, 0.8f, 0.1f);
	}

	private static final String DIR = "/home/dmyersturnbull/data/Bourne Lab/min_metric_change/";
	private static final File STATS_OUTPUT = new File("/home/dmyersturnbull/desktop/min_metric_change_data.txt");

	public static void run(float start, float stop, float step) throws FileNotFoundException {
		double max = 0;
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(STATS_OUTPUT);
			float argmax = -1;
			float value = start;
			Matrix top = null;
			while (value <= stop) {
				Sample sample = null;
				try {
					sample = runOn(value, new File(DIR + value + "_census.xml"), new File(DIR + value + "_sample.xml"));
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
				SimpleErrorMatrix matrix = new SimpleErrorMatrix();
				matrix.setOrderer(OrderDeterminationFactory.simpleWithScrew(2.0, 7.0));
				matrix.run(sample);
				double test = matrix.getDiagonalSum();
				pw.println("-------------------------" + value + "=" + matrix.getDiagonalSum() + "-------------------------");
				pw.println(matrix.getMatrix());
				pw.println("-----------------------------------------------------------------------");
				if (test > max) {
					max = test;
					argmax = value;
					top = matrix.getMatrix();
				}
				value += step;
			}
			pw.println("==============================TOP===================================");
			pw.println(argmax);
			pw.println(max);
			pw.println(top);
			pw.flush();
		} finally {
			if (pw != null) pw.close();
		}
	}

	public static Sample runOn(final float value, File censusFile, File sampleFile) throws IOException {
		NamesCensus census = new NamesCensus(0, NamesCensus.readNames(new File("src/main/resources/domain_symm_benchmark_names.list")));
		census.setAlgorithm(new AlgorithmGiver() {
			@Override
			public StructureAlignment getAlgorithm() {
				CeSymm ceSymm = new CeSymm();
				ceSymm.setMinimumMetricChange(value);
				return ceSymm;
			}
		});
		census.setOutputWriter(censusFile);
		census.run();
		census = null;
		SampleBuilder.buildSample(censusFile, sampleFile, new File("src/main/resources/domain_symm_benchmark.tsv"));
		return Sample.fromXML(sampleFile);
	}
}
