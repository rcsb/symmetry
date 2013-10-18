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
 * Static utilities related to math, particularly linear algebra.
 * @author dmyersturnbull
 */
public class MathUtils {

	public static RealMatrix readMatrix(File file) throws IOException {
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
		return matrix;
	}

	public static RealVector mapToShiftedVector(Map<Integer,Integer> map, int size) {
		RealVector vector = new ArrayRealVector(size);
		for (Map.Entry<Integer, Integer> entry : map.entrySet()) {
			vector.setEntry(entry.getKey()-1, entry.getValue());
		}
		return vector;
	}
	
	public static int argmax(double[] vector) {
		double max = 0;
		int argmax = -1;
		for (int i = 0; i < vector.length; i++) {
			double test = vector[i];
			if (test > max) {
				max = test;
				argmax = i;
			}
		}
		return argmax;
	}
	
	public static void printMatrix(RealMatrix matrix, int nDigits) throws IOException {
		StringBuilder sb = new StringBuilder();
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(nDigits);
		nf.setMinimumFractionDigits(nDigits);
		for (int i = 0; i < matrix.getRowDimension(); i++) {
			for (int j = 0; j < matrix.getColumnDimension(); j++) {
				sb.append(nf.format(matrix.getEntry(i, j)));
				if (j < matrix.getColumnDimension() - 1) sb.append("\t");
			}
			if (i < matrix.getRowDimension() - 1) sb.append(StatUtils.NEWLINE);
		}
		System.out.println(sb.toString());
	}

	public static int argmax(RealVector vector) {
		return argmax(vector.toArray());
	}

}
