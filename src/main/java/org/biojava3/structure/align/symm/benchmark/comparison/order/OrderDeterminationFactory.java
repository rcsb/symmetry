package org.biojava3.structure.align.symm.benchmark.comparison.order;

import org.biojava3.structure.align.symm.census2.Result;

/**
 * A factory for {@link OrderDetermination OrderDeterminations}.
 * @author dmyersturnbull
 */
public class OrderDeterminationFactory {

	public static OrderDetermination simpleWithScrew(final double tau, final double maxScrew) {
		return new OrderDetermination() {
			@Override
			public int getOrder(Result result) {
				if (result.getAxis() == null) return 1;
				int orderMethod1 = result.getOrder()==null||result.getOrder()==0? 1 : result.getOrder();
				int orderMethod2 = result.getAxis().guessOrder(tau, 8);
				double screw = result.getAxis().getScrew();
				if (screw > maxScrew) return 1;
				if (orderMethod2 != 1) return orderMethod2;
				return orderMethod1;
			}
		};
	}

	/**
	 * 
	 * @param theta
	 * @param threshold
	 * @return
	 */
	@Deprecated
	public int guessOrderFromAngle(double theta, double threshold) {
		final int maxOrder = 8;
		double bestDelta = threshold;
		int bestOrder = 1;
		for (int order = 2; order < maxOrder; order++) {
			double delta = Math.abs(2 * Math.PI / order - theta);
			System.out.println(delta);
			if (delta < bestDelta) {
				bestOrder = order;
				bestDelta = delta;
			}
		}
		return bestOrder;
	}
}
