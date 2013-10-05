package org.biojava3.structure.align.symm.benchmark.comparison.order;

import java.io.File;
import java.io.IOException;

import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

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
			System.out.println(StatUtils.formatD(i) + "\t" + StatUtils.formatD(test));
			i += step;
		}
		System.out.println("===========================================================");
		System.out.println(StatUtils.formatD(argmax) + "\t" + StatUtils.formatD(max));
		System.out.println(top);
	}
	
}
