package org.biojava.nbio.structure.align.symm.benchmark.comparison.order;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.benchmark.Case;
import org.biojava.nbio.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;
import org.biojava3.structure.align.symm.census3.CensusSymmetryGroup;
import org.biojava3.structure.align.symm.census3.stats.CensusStatUtils;
import org.biojava3.structure.align.symm.order.OrderDetector;
import org.biojava3.structure.align.symm.order.PeakCountingOrderDetector;
import org.biojava3.structure.align.symm.order.RotationOrderDetector;
import org.biojava3.structure.align.symm.order.SequenceFunctionOrderDetector;

/**
 * Test the performance of a {@link OrderDetector} against a benchmark.
 * @author dmyersturnbull
 */
public class BenchmarkOrderDetector {

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length >2) {
			System.err.println("Usage: " + BenchmarkOrderDetector.class.getSimpleName() + " input-benchmark-file [output-file]");
			return;
		}
		String inputFilename = args[0];
		String outputFilename = null;
		if(args.length>1) {
			outputFilename = args[1];
		}


		OrderDetector detector;
		detector = new RotationOrderDetector();
		//detector = new AngleOrderDetector(0.1);
		detector = new SequenceFunctionOrderDetector();

		// Read in benchmark
		AtomCache cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);
		ScopFactory.setScopDatabase(ScopFactory.VERSION_2_0_1, true);
		Sample sample = Sample.fromXML(new File(inputFilename));
		BenchmarkOrderDetector b = new BenchmarkOrderDetector(cache);

		//Set up output
		PrintStream out = System.out;
		if(outputFilename != null) {
			out = new PrintStream(outputFilename);
	}
	
		//Recalculate orders
		b.calculateOrder(sample, detector, out);

		//Calculate accuracy
		b.calculateOrder(sample, new PeakCountingOrderDetector(8),out);
		double acc = b.getOrderAccuracy(sample);
		System.out.println("Accuracy\t"+CensusStatUtils.formatP(acc));

		// Finish up
		if(out != System.out) {
			out.close();
		}
	}

	private AtomCache cache;

	public BenchmarkOrderDetector(AtomCache cache) {
		this.cache = cache;
	}

	/**
	 * Recalculates the order for each result in sample based on the specified
	 * OrderDetector
	 * @param sample
	 * @param detector
	 * @param out An output stream to report results (may be null)
	 */
	public void calculateOrder(Sample sample, OrderDetector detector, PrintStream out) {
		OrderDetermination determination = new OrderDetectorDeterminationAdaptor(detector, cache);
		for (Case c : sample.getData()) {
			CensusResult result = c.getResult();
			int order = determination.getOrder(result);
			result.setGroup(new CensusSymmetryGroup("C" + order));
			if(out != null) {
				out.println(result.getId()+"\t"+order);
			}
		}
	}

	/**
	 * Calculates the accuracy of the order detection
	 * @param sample
	 * @return
	 */
	public double getOrderAccuracy(Sample sample ) {
		OrderAccuracy acc = new OrderAccuracy(sample, CensusSignificanceFactory.forCeSymmTm(), GroupComparisonFactory.exact());
		return acc.getAccuracy();
	}

	/**
	 * Get the accuracy of a specified OrderDetector against a sample
	 * @param sample the benchmark, annotated with true results
	 * @param detector the order detection method
	 * @return the accuracy (TP+TN)/N
	 */
	public double testOrderDetector(Sample sample, OrderDetector detector) {
		calculateOrder(sample,detector,null);
		return getOrderAccuracy(sample);
	}
}
