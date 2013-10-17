package org.biojava3.structure.align.symm.order;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.RotationAxis;

/**
 * A more intelligent order-detection that uses angle, screw vector magnitude, and Spencer's method.
 * @author dmyersturnbull
 */
public class MultiMethodOrderDetector implements OrderDetector {

	private final double maxScrew;
	private final double angleError;

	public MultiMethodOrderDetector(double maxScrew, double angleError) {
		super();
		this.maxScrew = maxScrew;
		this.angleError = angleError;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {
		try {
			RotationAxis axis = new RotationAxis(afpChain);
			OrderDetector method1 = new SequenceFunctionOrderDetector();
			int orderMethod1 = method1.calculateOrder(afpChain, ca);
			OrderDetector method2 = new AngleOrderDetector(angleError);
			int orderMethod2 = method2.calculateOrder(afpChain, ca);
			double screw = (float) (Calc.amount(axis.getScrewTranslation()) / Calc.amount(axis.getRotationAxis()));
			if (screw > maxScrew) return 1;
			if (orderMethod2 != 1) return orderMethod2;
			return orderMethod1;
		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
	}

}
