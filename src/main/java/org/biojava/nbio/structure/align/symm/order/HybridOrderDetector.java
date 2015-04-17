package org.biojava.nbio.structure.align.symm.order;

import static java.lang.Math.*;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.RotationAxis;

/**
 * Hybrid method using AngleOrderDetectorPlus for the initial detection,
 * followed by checking possible higher orders via superposition
 * @author Spencer Bliven
 *
 */
public class HybridOrderDetector implements OrderDetector {
	
	// list of small primes
	private static final int[] primes = new int[] {2,3,5,7,11,13,17,19,23,29,31};

	private int maxOrder;
	private double error;
	private boolean normalizeError;
	
	public HybridOrderDetector(int maxOrder,double error) {
		this.maxOrder = maxOrder;
		this.error = error;
		this.normalizeError = true;

		if(maxOrder > primes[primes.length-1]) {
			throw new IllegalArgumentException("Orders greater than "+primes[primes.length-1]+ " are not supported.");
		}
	}
	
	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca)
			throws OrderDetectionFailedException {
		List<Integer> compatible = compatibleOrders(afpChain, ca);
		
		
		
		return -1;//TODO STUB
	}

	private List<Integer> compatibleOrders(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {
		// order -> probability
		List<Integer> compatible = new ArrayList<Integer>();
		RotationAxis axis;
		try {
			axis = new RotationAxis(afpChain);
		} catch (StructureException e) {
			throw new OrderDetectionFailedException(e);
		}
		double theta = axis.getAngle();

		for (int order = 1; order <= maxOrder; order++) {
			// Triangle wave starting at 0 with period 2pi/order
			double delta = abs(abs(theta*order/(2*PI)-.5)%1.0 - .5);
			// Triangle waves have amplitude 1, so need to un-normalize
			if(!normalizeError)
				delta *= 2*PI/order;

			if( delta <= error ) {
				compatible.add(order);
			}
		}
		return compatible;
	}

}
