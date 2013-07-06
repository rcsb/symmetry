package org.biojava3.structure.align.symm.benchmark.comparison;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava3.structure.align.symm.benchmark.Case;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

/**
 * A class to determine the accuracy of CE-Symm for determining order of rotational symmetry.
 * @author dmyerstu
 * @see AccuracyFinder, which considers only the binary choice of symmetric/asymmetric
 */
public class OrderAccuracy {

	private static final Logger logger = LogManager.getLogger(OrderAccuracy.class.getName());

	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + OrderAccuracy.class.getSimpleName() + " input-file");
			return;
		}
		File input = new File(args[0]);
		OrderAccuracy finder = new OrderAccuracy(input, SignificanceFactory.rotationallySymmetric(), GroupGuesserFactory.withMultiplesOk());
		System.out.println(finder);
	}

	public OrderAccuracy(File input, Significance sig, GroupGuesser guesser) throws IOException {
		this(Sample.fromXML(input), sig, guesser);
	}

	private int correct = 0;
	private int total = 0;

	public int guessOrderFromAngle(double theta, double threshold) {
		final int maxOrder = 8;
		double bestDelta = threshold;
		int bestOrder = 1;
		for (int order = 2; order < maxOrder; order++) {
			double delta = Math.abs(2 * Math.PI / order - theta);
			System.out.println(delta);
			if (delta < bestDelta) {
				bestOrder = order;
				bestDelta = delta;
			}
		}
		return bestOrder;
	}
	
	public OrderAccuracy(Sample sample, Significance sig, GroupGuesser guesser) {
		for (Case c : sample.getData()) {
			try {
				if (!c.getKnownInfo().hasRotationalSymmetry()) continue;
				if (!sig.isSignificant(c.getResult())) continue;
				Integer known = c.getKnownOrder();
				if (known == null) known = 1;
				Integer guess = c.getOrder();
				if (guess == null) guess = 1;
				boolean equiv = guesser.hasEquivalentOrder(known, guess);
				if (equiv) correct++;
				total++;
			} catch (RuntimeException e) {
				logger.error("Encountered an error on " + c.getScopId(), e);
			}
		}
	}

	@Override
	public String toString() {
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(2);
		return nf.format((double) correct / (double) (total) * 100.0) + "%";
	}

}
