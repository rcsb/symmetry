package org.biojava3.structure.align.symm.census2.stats.order;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

/**
 * A {@link CensusDecider} that treats any <em>divisor</em> of an order as potentially being the order itself.
 */
public class ErrorKernelDecider implements ConsensusDecider {

	private RealMatrix kernel;

	public ErrorKernelDecider(RealMatrix kernel) {
		this.kernel = kernel;
	}

	@Override
	public int decide(Map<Integer,Integer> countsByOrders) {
		RealVector counts = new ArrayRealVector(8);
		for (Map.Entry<Integer, Integer> entry : countsByOrders.entrySet()) {
			counts.setEntry(entry.getKey()-1, entry.getValue());
		}
		return decide(counts);
	}
	
	public int decide(RealVector counts) {

		double[] flows = new double[] {0, 0, 0, 0, 0, 0, 0, 0};

		for (int i = 1; i <= 8; i++) { // correct
			for (int j = 1; j <= 8; j++) { // putative TODO <=

				// the number of states we're allowed to use
				int m = Math.max(i, j) / Math.min(i, j) - 1;

				/*
				 *  use Chapman–Kolmogorov equation to find m-state transition kernel
				 *  We want the m-state transition because we're only concerned with flow from i to j,
				 *  which we know requires EXACTLY m steps
				 */
				RealMatrix mStateTransition = MatrixUtils.createRealIdentityMatrix(8);
				for (int k = 0; k < m; k++) mStateTransition = mStateTransition.multiply(kernel);

				// get the number of "j"s, our putative actual "i"s
				double count = counts.getEntry(j-1);

				/*
				 * We want to make "putative" flow into "correct"
				 * Ex: Matrix[4,2] = 0.1, the probability we got 2 instead of 4
				 * Then we're using mStateTransition[4, 2] = 0.1
				 * And thus we allow 0.1 * count to flow from 2 into 4
				 * So this is the correct indexing
				 */
				double flow = mStateTransition.getEntry(i-1, j-1); // TODO should be mStateTransition
				flows[i-1] += count * flow;

			}
		}

		double maxFlow = 0;
		int maximizingOrder = 0;
		for (int i = 1; i <= 8; i++) {
			double flow = flows[i-1];
			if (flow > maxFlow) {
				maxFlow = flow;
				maximizingOrder = i;
			}
		}
		return maximizingOrder;
	}

	public static ErrorKernelDecider fromMatrixFile() {
		try {
			return fromVectorFile(new File("src/main/resources/error_kernel.matrix"));
		} catch (IOException e) {
			throw new RuntimeException("Couldn't load matrix", e);
		}
	}

	public static ErrorKernelDecider fromVectorFile(File file) throws IOException {
		BufferedReader br = null;
		RealMatrix matrix = MatrixUtils.createRealMatrix(8, 8);
		try {
			br = new BufferedReader(new FileReader(file));
			String line = "";
			int r = 0;
			while ((line = br.readLine()) != null) {
				String[] parts = line.split("\t");
				RealVector vector = new ArrayRealVector(8);
				for (int i = 0; i < 8; i++) {
					vector.setEntry(i, Double.parseDouble(parts[i]));
				}
				matrix.setRow(r, vector.toArray());
				r++;
			}
		} finally {
			if (br != null) br.close();
		}
		printMatrix(matrix);
		return new ErrorKernelDecider(matrix);
	}

	private static void printMatrix(RealMatrix matrix) throws IOException {
		StringBuilder sb = new StringBuilder();
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(6);
		for (int i = 0; i < matrix.getRowDimension(); i++) {
			for (int j = 0; j < matrix.getColumnDimension(); j++) {
				sb.append(nf.format(matrix.getEntry(i, j)));
				if (j < matrix.getColumnDimension() - 1) sb.append("\t");
			}
			if (i < matrix.getRowDimension() - 1) sb.append(StatUtils.NEWLINE);
		}
		System.out.println(sb.toString());
	}

}
