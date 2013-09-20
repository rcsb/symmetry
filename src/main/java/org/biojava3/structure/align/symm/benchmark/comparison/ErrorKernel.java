package org.biojava3.structure.align.symm.benchmark.comparison;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.primes.Primes;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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

	private static final Logger logger = LogManager.getLogger(ErrorKernel.class.getName());

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2) {
			System.err.println("Usage: " + ErrorKernel.class.getSimpleName()
					+ " input-sample-file [matrix-output-file]");
			return;
		}
		ErrorKernel kernel = new ErrorKernel();
		Sample sample = Sample.fromXML(new File(args[0]));
		kernel.build(sample);
		System.out.println(kernel);
		System.out.println(kernel.vectorToString());
		Eigenpair steadyState = new Eigenpair(0, kernel.calcSteadyState(0.0000000001));
		System.out.println(steadyState);
		if (args.length > 1) {
			kernel.print(new File(args[1]));
		}

	}

	private double epsilon = 0.0000000001;

	public void setEpsilon(double epsilon) {
		this.epsilon = epsilon;
	}

	private RealMatrix kernel;
	private Significance significance = SignificanceFactory.forCeSymmTm();

	private double[] singleStepInverseMistakeRates = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	/**
	 * This is our function that maps divisors to probabilities. It is indexed
	 * starting at divisor=1. The transition of a composite number greater than
	 * 1 is defined to be zero (otherwise we'd be including multi-state
	 * transition probabilities).
	 */
	private double[] singleStepMistakeRates = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	public void build(Sample sample) {

		/*
		 * We'll use these to normalize the vectors
		 */
		int[] singleStepInverseMistakesPossible = new int[] { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		int[] singleStepMistakesPossible = new int[] { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

		Map<Integer, Integer> knownOrderCounts = new HashMap<Integer, Integer>();
		for (Case c : sample.getData()) {

			if (!significance.isSignificant(c.getResult())) continue;

			Integer order = c.getOrder();
			if (order == null || order < 1) order = 1;
			int knownOrder = c.getKnownOrder();
			StatUtils.plus(knownOrderCounts, knownOrder);

			/*
			 * Record all the possible primes, which are given by:
			 * phi(knownOrder / 1) for forward errors phi(order / 1) for inverse
			 * errors
			 */
			if (knownOrder > 1) {
				List<Integer> primeFactors = Primes.primeFactors(knownOrder);
				for (int primeFactor : primeFactors) {
					try {
						singleStepMistakesPossible[primeFactor - 1]++;
					} catch (ArrayIndexOutOfBoundsException e) {
						throw new RuntimeException("Got prime factor " + primeFactor
								+ ", which is too big for our vector");
					}
				}
			}
			if (order > 1) {
				List<Integer> primeFactors = Primes.primeFactors(order);
				for (int primeFactor : primeFactors) {
					try {
						singleStepInverseMistakesPossible[primeFactor - 1]++;
					} catch (ArrayIndexOutOfBoundsException e) {
						throw new RuntimeException("Got prime factor " + primeFactor
								+ ", which is too big for our vector");
					}
				}
			}

			/*
			 * If we get 1 instead of 12, we know that we made exactly 3
			 * mistakes: 2, 2, and 3 Therefore, we want to add 2 mistakes to
			 * divisor=2 and 1 to divisor=3
			 */
			int divisor = 0;
			double[] vector = null;
			if (knownOrder % order == 0) {
				divisor = knownOrder / order;
				vector = singleStepMistakeRates;
			} else if (order % knownOrder == 0) {
				divisor = order / knownOrder;
				vector = singleStepInverseMistakeRates;
			}
			if (divisor > 1) { // primeFactors fails when given 1
				List<Integer> primeFactors = Primes.primeFactors(divisor);
				for (int primeFactor : primeFactors) {
					try {
						vector[primeFactor - 1]++;
					} catch (ArrayIndexOutOfBoundsException e) {
						throw new RuntimeException("Got prime factor " + primeFactor
								+ ", which is too big for our vector");
					}
				}
			}
		}

		/*
		 * Normalize the mistake vectors
		 */
		for (int i = 1; i < singleStepMistakeRates.length; i++) {
			if (!Primes.isPrime(i + 1)) continue;
			if (singleStepMistakesPossible[i] == 0) throw new RuntimeException("Bad " + i);
			singleStepMistakeRates[i] = singleStepMistakeRates[i] / singleStepMistakesPossible[i];
		}
		for (int i = 1; i < singleStepInverseMistakeRates.length; i++) {
			if (!Primes.isPrime(i + 1)) continue;
			if (singleStepInverseMistakesPossible[i] == 0) throw new RuntimeException("Bad " + i);
			singleStepInverseMistakeRates[i] = singleStepInverseMistakeRates[i] / singleStepInverseMistakesPossible[i];
		}

		/*
		 * Now build the matrix from the vector This should be straightforward
		 */
		kernel = MatrixUtils.createRealMatrix(8, 8); // zeroes matrix;

		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				int knownOrder = i + 1;
				int order = j + 1;
				int divisor = 0;
				double[] vector = null;
				if (knownOrder % order == 0) {
					divisor = knownOrder / order;
					vector = singleStepMistakeRates;
				} else if (order % knownOrder == 0) {
					divisor = order / knownOrder;
					vector = singleStepInverseMistakeRates;
				}
				if (divisor != 0) {
					double flow = 0;
					if (divisor < vector.length - 1) {
						flow = vector[divisor - 1];
					} else if (Primes.isPrime(divisor)) { // we should always
						// include primes
						logger.error("Divisior " + divisor + " out of range");
					} // else it's just a composite number
					kernel.setEntry(i, j, flow + epsilon);
				} else {
					kernel.setEntry(i, j, epsilon);
				}
			}
		}


		/*
		 * Set the diagonal entries as 1 minus the rest of the row
		 */
		for (int i = 1; i <= 8; i++) {
			double n = 1;
			for (int j = 1; j <= 8; j++) {
				n -= kernel.getEntry(i - 1, j - 1);
			}
			kernel.setEntry(i - 1, i - 1, n);
		}

	}

	public static class Eigenpair {
		private double eigenvalue;
		private double[] eigenvector;
		public double getEigenvalue() {
			return eigenvalue;
		}
		public void setEigenvalue(double eigenvalue) {
			this.eigenvalue = eigenvalue;
		}
		public double[] getEigenvector() {
			return eigenvector;
		}
		public void setEigenvector(double[] eigenvector) {
			this.eigenvector = eigenvector;
		}
		public Eigenpair(double eigenvalue, double[] eigenvector) {
			this.eigenvalue = eigenvalue;
			this.eigenvector = eigenvector;
		}

		public Eigenpair(double eigenvalue, RealVector eigenvector) {
			this(eigenvalue, eigenvector.toArray());
		}
		@Override
		public String toString() {
			NumberFormat nf = new DecimalFormat();
			nf.setMaximumFractionDigits(5);
			if (eigenvector == null) return "null";
			StringBuilder sb = new StringBuilder();
			sb.append(StatUtils.formatD(eigenvalue));
			sb.append(": [");
			for (int i = 0; i < eigenvector.length; i++) {
				sb.append(nf.format(eigenvector[i]));
				if (i < eigenvector.length - 1) sb.append(", ");
			}
			sb.append("]");
			return sb.toString();
		}
	}

	public double[] calcSteadyState(double minChange) {
		RealMatrix m = kernel;
		double delta = Double.POSITIVE_INFINITY;
		while (delta > minChange) {
			double prev = m.getFrobeniusNorm();
			m = m.multiply(kernel);
			double next = m.getFrobeniusNorm();
			if (Math.abs(next - prev) < minChange) {
				break;
			}
		}
		return m.getRow(0); // arbitrary row
	}

	public double[] calcSteadyState(int power) {
		RealMatrix m = kernel.power(power);
		return m.getRow(0); // arbitrary row
	}

	/**
	 * Calculates a steady-state distribution by eigendecomposition.
	 * @return An Eigenpair containing the stationary distribution and the the dominant eigenvalue, which should be equal to 1.
	 * TODO Why doesn't this work?
	 */
	public Eigenpair calcDominantEigenpair() {
		EigenDecomposition decomposition = new EigenDecomposition(kernel);
		RealMatrix eigenvalues = decomposition.getD();
		RealVector pi = null;
		double lambda = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < kernel.getRowDimension(); i++) {
			if (Math.abs(eigenvalues.getEntry(i, i)) > lambda) {
				lambda = eigenvalues.getEntry(i, i);
				pi = decomposition.getEigenvector(i);
			}
		}
		return new Eigenpair(lambda, pi);
	}

	public RealMatrix getKernel() {
		return kernel;
	}

	public double[] getVector() {
		return singleStepMistakeRates;
	}

	public void print(File output) throws IOException {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new FileWriter(output));
			print(pw);
		} finally {
			if (pw != null) pw.close();
		}
	}

	public void print(PrintWriter output) {
		//		kernel.print(output, 5, 4);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				sb.append((kernel.getEntry(i, j)));
				if (j < 7) sb.append(" ");
			}
			sb.append(StatUtils.NEWLINE);
		}
		return sb.toString();
	}

	public String vectorToString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for (int i = 0; i < singleStepMistakeRates.length; i++) {
			sb.append(StatUtils.formatD(singleStepMistakeRates[i]));
			if (i < singleStepMistakeRates.length - 1) sb.append(" ");
		}
		sb.append("]" + StatUtils.NEWLINE);
		sb.append("[");
		for (int i = 0; i < singleStepInverseMistakeRates.length; i++) {
			sb.append(StatUtils.formatD(singleStepInverseMistakeRates[i]));
			if (i < singleStepInverseMistakeRates.length - 1) sb.append(" ");
		}
		sb.append("]" + StatUtils.NEWLINE);
		return sb.toString();
	}

}
