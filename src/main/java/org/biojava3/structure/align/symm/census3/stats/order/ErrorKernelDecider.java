package org.biojava3.structure.align.symm.census3.stats.order;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * A {@link CensusDecider} that treats any <em>divisor</em> of an order as potentially being the order itself.
 * @author dmyersturnbull
 * @deprecated This was being used, but now we're focusing on trying to improve the order-detector instead
 */
public class ErrorKernelDecider implements ConsensusDecider {

	private RealMatrix kernel;

	public ErrorKernelDecider(RealMatrix kernel) {
		this.kernel = kernel;
	}

	@Override
	public int decide(Map<Integer,Integer> countsByOrders) {
		return decide(MathUtils.mapToShiftedVector(countsByOrders, 8));
	}
	
	public int decide(RealVector counts) {

		RealMatrix upper = MatrixUtils.createRealMatrix(8, 8);
		RealMatrix lower = MatrixUtils.createRealMatrix(8, 8);
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j <= i; j++) {
				lower.setEntry(i, j, kernel.getEntry(i, j));
			}
			for (int j = i + 1; j < 8; j++) {
				upper.setEntry(i, j, kernel.getEntry(i, j));
			}
		}
		
		double[] flows = new double[] {0, 0, 0, 0, 0, 0, 0, 0};

		for (int i = 1; i <= 8; i++) { // correct
			for (int j = 1; j <= 8; j++) { // putative

				// the number of states we're allowed to use
				int m = Math.max(i, j) / Math.min(i, j) - 1;
				
				final RealMatrix half = i>j? upper : lower;
				final double count = i>j? -counts.getEntry(i-1) : counts.getEntry(j-1); // note the - sign

				/*
				 *  use Chapman–Kolmogorov equation to find m-state transition kernel
				 *  We want the m-state transition because we're only concerned with flow from i to j,
				 *  which we know requires EXACTLY m steps
				 */
				RealMatrix mStateTransition = MatrixUtils.createRealIdentityMatrix(8);
				for (int k = 0; k < m; k++) mStateTransition = mStateTransition.multiply(half);

				/*
				 * We want to make "putative" flow into "correct"
				 * Ex: Matrix[4,2] = 0.1, the probability we got 2 instead of 4
				 * Then we're using mStateTransition[4, 2] = 0.1
				 * And thus we allow 0.1 * count to flow from 2 into 4
				 * So this is the correct indexing
				 */
				double flow = mStateTransition.getEntry(i-1, j-1);
				flows[i-1] += count * flow;

			}
		}

		return MathUtils.argmax(flows) + 1;
	}

	public static ErrorKernelDecider fromMatrixFile() {
		try {
			return fromMatrixFile(new File("src/main/resources/error_kernel.matrix"));
		} catch (IOException e) {
			throw new RuntimeException("Couldn't load matrix", e);
		}
	}

	public static ErrorKernelDecider fromMatrixFile(File file) throws IOException {
		RealMatrix matrix = MathUtils.readMatrix(file);
		MathUtils.printMatrix(matrix, 6);
		return new ErrorKernelDecider(matrix);
	}

}
