package org.biojava3.structure.align.symm.benchmark.comparison.order;

import java.io.File;
import java.io.IOException;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.align.symm.benchmark.Case;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;
import org.biojava3.structure.align.symm.order.OrderDetector;
import org.biojava3.structure.align.symm.order.RotationOrderDetector;

/**
 * Test the performance of a {@link OrderDetector} against a benchmark.
 * @author dmyersturnbull
 */
public class BenchmarkOrderDetector {

	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + BenchmarkOrderDetector.class.getSimpleName() + " input-benchmark-file");
			return;
		}
		Sample sample = Sample.fromXML(new File(args[0]));
		BenchmarkOrderDetector b = new BenchmarkOrderDetector(new AtomCache());
		double acc = b.testOrderDetector(sample, new RotationOrderDetector());
		System.out.println("accuracy=" + StatUtils.formatP(acc));
	}
	
	private AtomCache cache;

	public BenchmarkOrderDetector(AtomCache cache) {
		this.cache = cache;
	}

	public double testOrderDetector(Sample sample, OrderDetector detector) {
		OrderDetermination determination = new OrderDetectorDeterminationAdaptor(detector, cache);
		for (Case c : sample.getData()) {
			Result result = c.getResult();
			result.setOrder(determination.getOrder(result));
		}
		OrderAccuracy acc = new OrderAccuracy(sample, SignificanceFactory.generallySymmetric(), GroupComparisonFactory.exact());
		return acc.getAccuracy();
	}

	
}
