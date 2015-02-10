package org.biojava.nbio.structure.align.symm.benchmark.comparison.order;

import java.io.File;
import java.io.IOException;

import org.biojava.nbio.structure.align.symm.benchmark.Sample;
import org.biojava.nbio.structure.align.symm.census3.stats.CensusStatUtils;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Test the <em>performance</em> of order-detection methods.
 * @author dmyersturnbull
 * TODO Improve for new functionality/tests
 */
public class OrderDeterminationTester {

	public static void main(String[] args) throws IOException {
		run(Sample.fromXML(new File(args[0])), 0.0, 20, 0.001);
	}

	// tau = 6.7 is top
	// angle = 1.91°
	private static final double BEST_TAU = 6.7;
	
	public static void run(Sample sample, double start, double stop, double step) {
		double argmax = -1;
		double max = 0;
		Matrix top = null;
		double i = start;
		while (i <= stop) {
			SimpleErrorMatrix mx = new SimpleErrorMatrix();
			mx.setOrderer(OrderDeterminationFactory.simpleWithScrew(i * Math.PI / 180, BEST_TAU));
			mx.run(sample);
			double test = mx.getDiagonalSum();
			if (test > max) {
				max = test;
				argmax = i;
				top = mx.getMatrix();
			}
			System.out.println(CensusStatUtils.formatD(i) + "\t" + CensusStatUtils.formatD(test));
			i += step;
		}
		System.out.println("===========================================================");
		System.out.println(CensusStatUtils.formatD(argmax) + "\t" + CensusStatUtils.formatD(max));
		System.out.println(top);
	}
	
}
