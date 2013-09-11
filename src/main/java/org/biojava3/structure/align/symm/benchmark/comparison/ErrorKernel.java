package org.biojava3.structure.align.symm.benchmark.comparison;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.primes.Primes;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.align.symm.benchmark.Case;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;
import org.biojava3.structure.utils.SymmetryTools;

/**
 * Uses a benchmark {@link Sample} to generate a transition kernel (as a matrix)
 * that gives the "1-step" probability that CE-Symm's order-detection (main
 * method; see {@link SymmetryTools}) gives a different answer. Here, "1-step"
 * means something very specific and is described below. We model the error as a
 * Markov chain in which each transition not on the diagonal of the kernel is a
 * decision by CE-Symm to use a divisor of the order. This makes a number of
 * assumptions:
 * <ol>
 * <li>CE-Symm only returns divisors of the correct order. This is certainly not
 * completely true, but it's close.</li>
 * <li>The error rate is a function only of the correct order divided by the
 * order. Therefore, the transition from C6 to C3 equals the transition from C8
 * to C4. This assumption is incorrect, but reasons that cause differences in
 * these transitions in real data are expected to be a result of real biological
 * phenomena or artifacts of the database (SCOP) rather than the mechanics of
 * CE-Symm or its order-detection.</li>
 * </ol>
 * The objective here is to get decent transition probabilities for cases where
 * we have few data points in the benchmark set; namely, we have too few points
 * for C8 to C4, for C8 to C2, for C6 to C3, for C6 to C2, and for C4 to C2.
 * Making the above assumptions allows us to determine these rates with less
 * random error.
 * 
 * @author dmyerstu
 */
public class ErrorKernel {

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2) {
			System.err.println("Usage: " + ErrorKernel.class.getSimpleName()
					+ " input-sample-file [matrix-output-file]");
			return;
		}
		Sample sample = Sample.fromXML(new File(args[0]));
		ErrorKernel kernel = new ErrorKernel(sample);
		System.out.println(kernel);
		if (args.length > 1) {
			kernel.print(new File(args[1]));
		}

	}

	private Matrix kernel;

	private Significance significance = SignificanceFactory.rotationallySymmetricSmart();

	public ErrorKernel(Sample sample) {

		/**
		 * This is our function that maps divisors to probabilities. It is
		 * indexed starting at divisor=2. The transition of a composite number greater than 1
		 * is defined to be zero (otherwise we'd be including multi-state transition
		 * probabilities).
		 */
		int[] singleStepMistakeRates = new int[] { 0, 0, 0, 0, 0 };

		Map<Integer, Integer> knownOrderCounts = new HashMap<Integer, Integer>();
		for (Case c : sample.getData()) {

			if (c.getKnownOrder() < 2)
				continue;
			if (c.getOrder() == null || c.getOrder() < 2)
				continue;
			if (!significance.isSignificant(c.getResult()))
				continue;

			int order = c.getOrder();
			int knownOrder = c.getKnownOrder();
			StatUtils.plus(knownOrderCounts, knownOrder);

			if (knownOrder % order == 0) {

				/*
				 * If we get 2 instead of 12, we know that we made exactly 3 mistakes:
				 * 2, 2, and 3
				 * Therefore, we want to add 2 mistakes to divisor=2 and 1 to divisor=3
				 */
				int divisor = knownOrder / order;
				if (divisor > 1) { // primeFactors fails when given 1
					List<Integer> primeFactors = Primes.primeFactors(divisor);
					for (int primeFactor : primeFactors) {
						singleStepMistakeRates[primeFactor]++;
					}
				} else {
					singleStepMistakeRates[1]++;
				}

			} else if (order % knownOrder == 0) {
				// okay, this case is weird, but we need to take it into account
				// TODO I'm not sure what to do with negative numbers though
			}

		}

		kernel = new Matrix(7, 7); // zeroes matrix;
		for (int i = 0; i < 7; i++) {
			double rowSum = 0;
			for (int j = 0; j < 7; j++) {
				int knownOrder = i + 2;
				int order = j + 2;
				if (knownOrder % order == 0) {
					int divisor = knownOrder / order;
					double flow = singleStepMistakeRates[divisor];
					kernel.set(i, j, flow);
					rowSum += flow;
				}
			}
			for (int j = 0; j < 7; j++) {
				kernel.set(i, j, kernel.get(i, j) / rowSum);
			}
		}

	}

	public Matrix getKernel() {
		return kernel;
	}

	public void print(File output) throws IOException {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new FileWriter(output));
			print(pw);
		} finally {
			if (pw != null)
				pw.close();
		}
	}

	public void print(PrintWriter output) {
		kernel.print(output, 5, 4);
	}

	@Override
	public String toString() {
		return kernel.toString();
	}

}
