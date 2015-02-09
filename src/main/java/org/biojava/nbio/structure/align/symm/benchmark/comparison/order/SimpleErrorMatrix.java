package org.biojava.nbio.structure.align.symm.benchmark.comparison.order;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.nbio.structure.align.symm.benchmark.Case;
import org.biojava.nbio.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census3.CensusSignificance;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;

/**
 * Builds a matrix indexed by the correct order and the order CE-Symm found.
 * Elements are the number of cases found.
 * 
 * @author dmyersturnbull
 */
public class SimpleErrorMatrix {

	private static final Logger logger = LogManager.getLogger(SimpleErrorMatrix.class.getName());

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2) {
			System.err.println("Usage: " + SimpleErrorMatrix.class.getSimpleName() + " input-sample-file");
			return;
		}
		Sample sample = Sample.fromXML(new File(args[0]));
		SimpleErrorMatrix matrix = new SimpleErrorMatrix();
		matrix.run(sample);
		System.out.println(matrix);
		System.out.println(matrix.getDiagonalSum());
	}

	private Matrix matrix = new Matrix(8, 8);

	private CensusSignificance significance = CensusSignificanceFactory.forCeSymmTm();
	private OrderDetermination orderer = OrderDeterminationFactory.simpleWithScrew(2 * Math.PI / 180, 7);

	public void setOrderer(OrderDetermination orderer) {
		this.orderer = orderer;
	}

	public void run(Sample sample) {
		for (Case c : sample.getData()) {
			if (!significance.isSignificant(c.getResult())) continue;
			int knownOrder = c.getKnownOrder();
			int order = orderer.getOrder(c.getResult());
			matrix.set(knownOrder - 1, order - 1, matrix.get(knownOrder - 1, order - 1) + 1);
			if (knownOrder % order != 0 && order % knownOrder != 0) {
//				logger.info("Found: " + knownOrder + "," + order);
			}
		}
	}

	public Matrix getMatrix() {
		return matrix;
	}

	public int getDiagonalSum() {
		int sum = 0;
		for (int i = 0; i < 8; i++) {
			sum += matrix.get(i, i);
		}
		return sum;
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
		matrix.print(output, 2, 0);
	}

	public void setSignificance(CensusSignificance significance) {
		this.significance = significance;
	}

	@Override
	public String toString() {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		print(pw);
		return sw.toString();
	}

}
