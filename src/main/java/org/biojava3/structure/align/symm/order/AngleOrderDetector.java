package org.biojava3.structure.align.symm.order;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.RotationAxis;

/**
 * Guesses an order of rotational symmetry from the angle.
 * @author dmyersturnbull
 */
public class AngleOrderDetector implements OrderDetector {

	private int maxOrder = 8;
	private final double angleError;

	public AngleOrderDetector(double angleError) {
		super();
		this.angleError = angleError;
	}

	public AngleOrderDetector(int maxOrder, double angleError) {
		super();
		this.maxOrder = maxOrder;
		this.angleError = angleError;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {
		
		try {
			RotationAxis axis = new RotationAxis(afpChain);
			double theta = axis.getAngle();

			double bestDelta = angleError;
			int bestOrder = 1;
			for (int order = 2; order < maxOrder; order++) {
				double delta = Math.abs(2 * Math.PI / order - theta);
				if (delta < bestDelta) {
					bestOrder = order;
					bestDelta = delta;
				}
			}
			return bestOrder;
		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
	}


}
