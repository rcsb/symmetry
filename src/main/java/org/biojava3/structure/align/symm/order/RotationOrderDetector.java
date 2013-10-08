package org.biojava3.structure.align.symm.order;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.jama.Matrix;

/**
 * Detects order by rotating the structure by angles that correspond to orders.
 * This one could be smart.
 * TODO Needs lots of work
 * @author dmyersturnbull
 */
public class RotationOrderDetector implements OrderDetector {

	private int maxOrder = 8;
	private final double minTmScore;

	public RotationOrderDetector(double minTmScore) {
		super();
		this.minTmScore = minTmScore;
	}

	public RotationOrderDetector(int maxOrder, double minTmScore) {
		super();
		this.maxOrder = maxOrder;
		this.minTmScore = minTmScore;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);

			AFPChain clone = (AFPChain) afpChain.clone();
			Atom[] ca2 = null;

			double bestTmScore = minTmScore;
			int argmax = 1;

			// TODO apply alignment
			
			for (int order = 1; order <= maxOrder; order++) {

				ca2 = StructureTools.cloneCAArray(ca); // reset rotation for new order
				Matrix rotation = axis.getRotationMatrix(2*Math.PI / order); // will apply repeatedly

				/*
				 * If C6, we should be able to rotate 6 times and still get a decent superposition
				 */
				double lowestTmScore = Double.POSITIVE_INFINITY;
				for (int j = 1; j < order; j++) {
					Calc.rotate(ca2, rotation); // rotate repeatedly
					double tmScore = AFPChainScorer.getTMScore(clone, ca, ca2);
					if (tmScore < lowestTmScore) {
						lowestTmScore = tmScore;
					}
				}

				if (lowestTmScore > bestTmScore) {
					bestTmScore = lowestTmScore;
					argmax = order;
				}

			}

			return argmax;

		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
	}

}
