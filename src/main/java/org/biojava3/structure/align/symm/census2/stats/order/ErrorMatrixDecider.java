package org.biojava3.structure.align.symm.census2.stats.order;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * A simplified version of {@link ErrorKernelDecider}.
 * @author dmyersturnbull
 */
public class ErrorMatrixDecider implements ConsensusDecider {

	private RealMatrix matrix;

	public ErrorMatrixDecider(RealMatrix matrix) {
		this.matrix = matrix;
	}

	@Override
	public int decide(Map<Integer,Integer> countsByOrders) {
		return decide(MathUtils.mapToShiftedVector(countsByOrders));
	}
	
	public int decide(RealVector counts) {
		RealVector modified = matrix.operate(counts);
		return MathUtils.argmax(modified) + 1;
	}

	public static ErrorMatrixDecider fromMatrixFile() {
		try {
			return fromMatrixFile(new File("src/main/resources/full_error_matrix.matrix"));
		} catch (IOException e) {
			throw new RuntimeException("Couldn't load matrix", e);
		}
	}

	public static ErrorMatrixDecider fromMatrixFile(File file) throws IOException {
		RealMatrix matrix = MathUtils.readMatrix(file);
		MathUtils.printMatrix(matrix, 6);
		return new ErrorMatrixDecider(matrix);
	}

}
