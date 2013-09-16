package org.biojava3.structure.align.symm.census2.stats.order;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import org.biojava.bio.structure.jama.Matrix;

/**
 * A {@link CensusDecider} that treats any <em>divisor</em> of an order as potentially being the order itself.
 */
public class ErrorKernelDecider implements ConsensusDecider {

	private Matrix kernel;

	public ErrorKernelDecider(Matrix kernel) {
		super();
		this.kernel = kernel;
	}

	@Override
	public int decide(Map<Integer,Integer> countsByOrders) {

		double[] flows = new double[] {0, 0, 0, 0, 0, 0, 0};

		for (int i = 2; i <= 8; i++) { // correct
			for (int j = 2; j <= 8; j++) { // putative

				// the number of states we're allowed to use
				int m = i / j - 1;

				/*
				 *  use Chapman–Kolmogorov equation to find m-state transition kernel
				 *  We want the m-state transition because we're only concerned with flow from i to j,
				 *  which we know requires EXACTLY m steps
				 */
				Matrix mStateTransition = Matrix.identity(7, 7);
				for (int k = 0; k < m; k++) mStateTransition = mStateTransition.times(kernel);

				// get the number of "j"s, our putative actual "i"s
				Integer count = countsByOrders.get(j);
				if (count == null) count = 0;

				/*
				 * We want to make "putative" flow into "correct"
				 * Ex: Matrix[4,2] = 0.1, the probability we got 2 instead of 4
				 * Then we're using mStateTransition[4, 2] = 0.1
				 * And thus we allow 0.1 * count to flow from 2 into 4
				 * So this is the correct indexing
				 */
				double flow = mStateTransition.get(i - 2, j - 2);
				flows[i-2] += count * flow;

			}
		}

		double maxFlow = 0;
		int maximizingOrder = 0;
		for (int i = 2; i <= 8; i++) {
			double flow = flows[i-2];
			if (flow > maxFlow) {
				maxFlow = flow;
				maximizingOrder = i;
			}
		}
		return maximizingOrder;
	}

	public static ErrorKernelDecider fromMatrixFile() {
		try {
			return fromMatrixFile(new File("src/main/resources/error_kernel.matrix"));
		} catch (IOException e) {
			throw new RuntimeException("Couldn't load matrix", e);
		}
	}

	public static ErrorKernelDecider fromMatrixFile(File file) throws IOException {
		Matrix kernel;
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(file));
			kernel = Matrix.read(br);
		} finally {
			if (br != null) br.close();
		}
		return new ErrorKernelDecider(kernel);
	}

}
