package org.biojava3.structure.align.symm.benchmark.comparison.order;

import java.io.File;
import java.io.IOException;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.benchmark.Case;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;
import org.biojava3.structure.align.symm.census3.CensusSymmetryGroup;
import org.biojava3.structure.align.symm.census3.stats.CensusStatUtils;
import org.biojava3.structure.align.symm.order.OrderDetector;
import org.biojava3.structure.align.symm.order.PeakCountingOrderDetector;
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
		ScopFactory.setScopDatabase(ScopFactory.VERSION_2_0_1);
		Sample sample = Sample.fromXML(new File(args[0]));
		AtomCache cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);
		BenchmarkOrderDetector b = new BenchmarkOrderDetector(cache);
		double acc = b.testOrderDetector(sample, new PeakCountingOrderDetector(8));
		System.out.println("accuracy=" + CensusStatUtils.formatP(acc));
	}
	
	private AtomCache cache;

	public BenchmarkOrderDetector(AtomCache cache) {
		this.cache = cache;
	}

	public double testOrderDetector(Sample sample, OrderDetector detector) {
		OrderDetermination determination = new OrderDetectorDeterminationAdaptor(detector, cache);
		for (Case c : sample.getData()) {
			CensusResult result = c.getResult();
			result.setGroup(new CensusSymmetryGroup("C" + determination.getOrder(result)));
		}
		OrderAccuracy acc = new OrderAccuracy(sample, CensusSignificanceFactory.forCeSymmTm(), GroupComparisonFactory.exact());
		return acc.getAccuracy();
	}

	
}
