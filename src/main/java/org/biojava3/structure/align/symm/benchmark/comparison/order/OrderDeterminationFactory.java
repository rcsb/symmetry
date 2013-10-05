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
	
}
