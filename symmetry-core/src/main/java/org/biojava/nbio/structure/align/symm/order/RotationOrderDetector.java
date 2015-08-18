package org.biojava.nbio.structure.align.symm.order;

import java.util.Arrays;

import org.apache.commons.math3.util.Pair;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.symmetry.internal.OrderDetector;
import org.biojava.nbio.structure.symmetry.internal.RefinerFailedException;

import static java.lang.Math.*;

/**
 * Detects order by analyzing the goodness of fit as the protein is rotated
 * around the axis of symmetry.
 * @author Spencer Bliven
 */
public class RotationOrderDetector implements OrderDetector {
	public static enum RotationOrderMethod {
		/**
		 * Model as a sum of sin^2 terms of decreasing period:
		 *   f(angle) = 0 + a1*sin(angle/2)^2 + a2*sin(2*angle/2)^2 + a3*sin(3*angle/2)^2 + ...
		 * No intercept, includes angles 0 to pi for fitting.
		 * Return the order with the highest amplitude.
		 */
		HARMONICS,
		/**
		 * Adds an intercept to the HARMONICS method:
		 *   f(angle) = a0 + a1*sin(angle/2)^2 + a2*sin(2*angle/2)^2 + a3*sin(3*angle/2)^2 + ...
		 * Fits angles from pi/maxOrder to pi.
		 * Return the order with the highest amplitude.
		 */
		HARMONICS_FLOATING,
		/**
		 * Fits a function for each order independently:
		 *   f(angle, order) = a0 + ai*sin(i*angle/2)^2
		 * Fits angles from pi/maxOrder to pi.
		 * Return the order with the highest amplitude.
		 */
		SINGLE_HARMONIC_AMP,
		/**
		 * Fits a function for each order independently:
		 *   f(angle, order) = a0 + ai*sin(i*angle/2)^2
		 * Fits angles from pi/maxOrder to pi.
		 * Return the order with the lowest sum squared error.
		 */
		SINGLE_HARMONIC_SSE,
		/**
		 * Fits a poorly conceived cusp function for each order:
		 *   f(angle, order) = a0 + ai*sqrt( 1 - cos(order*angle) );
		 * Fits angles from pi/maxOrder to pi.
		 * Return the order with the highest amplitude.
		 */
		SINGLE_CUSP_AMP,
		/**
		 * Fits a poorly conceived cusp function for each order:
		 *   f(angle, order) = a0 + ai*sqrt( 1 - cos(order*angle) );
		 * Fits angles from pi/maxOrder to pi.
		 * Return the order with the lowest sum squared error.
		 */
		SINGLE_CUSP_SSE,
		/**
		 * Fits a better derived cusp function for each order:
		 *   f(angle, order) = a0 + ai*sqrt( 2 - 2*cos(2*PI/order*triangle(angle,order)) )
		 * where triangle(angle,order) gives the triangle wave with period 2pi/order.
		 * Fits angles from pi/maxOrder to pi.
		 * Return the order with the highest amplitude.
		 */
		SINGLE_CUSP_FIXED_AMP,
		/**
		 * Fits a better derived cusp function for each order:
		 *   f(angle, order) = a0 + ai*sqrt( 2 - 2*cos(2*PI/order*triangle(angle,order)) )
		 * where triangle(angle,order) gives the triangle wave with period 2pi/order.
		 * Fits angles from pi/maxOrder to pi.
		 * Return the order with the lowest sum squared error.
		 */
		SINGLE_CUSP_FIXED_SSE,
	}

	private int maxOrder;

	public static final double DEFAULT_ANGLE_INCR = Math.toRadians(5);
	private double angleIncr = DEFAULT_ANGLE_INCR; // angular resolution

	private RotationOrderMethod method;
	private double minAngle;
	public RotationOrderDetector() {
		this(8);
	}

	public RotationOrderDetector(int maxOrder) {
		this(maxOrder,RotationOrderMethod.SINGLE_HARMONIC_SSE);
	}

	public RotationOrderDetector(int maxOrder, RotationOrderMethod method) {
		this.maxOrder = maxOrder;
		this.method= method;
		switch(method) {
		case HARMONICS:
			this.minAngle = 0;
			break;
		default:
			this.minAngle = PI/maxOrder;
		}
	}

	@Override
	public String toString() {
		return getClass().getSimpleName()+"[method="+method+",maxOrder="+maxOrder+"]";
	}
	public void setMethod(RotationOrderMethod m) {
		method = m;
	}
	public RotationOrderMethod getMethod() {
		return method;
	}
	public int getMaxOrder() {
		return maxOrder;
	}

	public void setMaxOrder(int maxOrder) {
		this.maxOrder = maxOrder;
	}

	public double getAngleIncr() {
		return angleIncr;
	}

	public void setAngleIncr(double angleIncr) {
		this.angleIncr = angleIncr;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws RefinerFailedException {
		//TODO only use aligned residues, rather than the whole ca
		try {

			RotationAxis axis = new RotationAxis(afpChain);

			// Use C1 order if the axis is undefined
			if(!axis.isDefined()) {
				return 1;
			}
			
			// Calculate optimum
			switch(method) {
			case HARMONICS: {
				double[] coefficients = tryAllOrders(ca, axis, false);
				
				return maxIndex(coefficients, 1);
			}
			case HARMONICS_FLOATING: {
				double[] coefficients = tryAllOrders(ca, axis, true);
				
				return maxIndex(coefficients, 1);
			}
			case SINGLE_HARMONIC_AMP:
			case SINGLE_CUSP_AMP:
			case SINGLE_CUSP_FIXED_AMP: {
				// Calculate weights for each order
				double[] coefficients = trySingleOrdersByAmp(ca,axis);

				// Find order with maximum weight
				return maxIndex(coefficients,0) + 1;
			}
			case SINGLE_HARMONIC_SSE:
			case SINGLE_CUSP_SSE:
			case SINGLE_CUSP_FIXED_SSE: {
				// Calculate weights for each order
				double[] coefficients = trySingleOrdersBySSE(ca,axis);
				
				// Find order with minimum SSE
				return minIndex(coefficients,0) + 1;
			}
			default:
				throw new UnsupportedOperationException("Unimplemented method "+method);
			}

		} catch (StructureException e) {
			throw new RefinerFailedException(e);
		}
	}

	/**
	 * Returns an array of {@link #superpositionDistance(Atom[], Atom[]) superposition distances} of rotations of {@code ca}.
	 * The {@code n}th element in the array corresponds to a rotation by {@code degreesIncrement * n} degrees.
	 */
	public static Pair<double[],double[]> sampleRotations(Atom[] ca, RotationAxis axis, double degreesIncrement) throws StructureException {
		final double angleIncr = Math.toRadians(degreesIncrement);
		final int steps = (int)floor(2*PI/angleIncr);

		double[] angles = new double[steps];
		double[] distances = new double[steps];

		Atom[] ca2 = StructureTools.cloneAtomArray(ca);
		double angle = 0;

		for (int step=0; step<steps;step++) {
			angles[step] = angle;
			double dist = superpositionDistance(ca, ca2);
			distances[step] = dist;
			// Rotate for next step
			axis.rotate(ca2, angleIncr);
			angle += angleIncr;
		}

		return new Pair<double[], double[]>(angles, distances);

	}

	/**
	 * Provide a rough alignment-free metric for the similarity between two
	 * superimposed structures.
	 *
	 * The average distance from each atom to the closest atom in the other
	 * is used.
	 * @param ca1 first structure
	 * @param ca2 second structure
	 * @return the average distance to the closest atom
	 * @throws StructureException if an error occurs finding distances between atoms
	 */
	public static double superpositionDistance(Atom[] ca1, Atom[] ca2) throws StructureException {

		// Store the closest distance yet found
		double[] bestDist1 = new double[ca1.length];
		double[] bestDist2 = new double[ca2.length];
		Arrays.fill(bestDist1, Double.POSITIVE_INFINITY);
		Arrays.fill(bestDist2, Double.POSITIVE_INFINITY);

		for(int i=0;i<ca1.length;i++) {
			for(int j=0;j<ca2.length;j++) {
				double dist = Calc.getDistanceFast(ca1[i], ca2[j]);
				if( dist < bestDist1[i]) {
					bestDist1[i] = dist;
				}
				if( dist < bestDist2[j]) {
					bestDist2[j] = dist;
				}
			}
		}

		double total = 0;
		for(int i=0;i<ca1.length;i++) {
			total += sqrt(bestDist1[i]);
		}
		for(int j=0;j<ca2.length;j++) {
			total += sqrt(bestDist2[j]);
		}

		double dist = total/(ca1.length+ca2.length);
		return dist;
	}

	protected double[] getAngles() {
		final double firstAngle = ceil(this.minAngle/angleIncr)*angleIncr; // first valid peak
		final double maxAngle = PI;
		// Number of angle steps
		final int steps = (int)floor((maxAngle-firstAngle)/angleIncr);

		double[] angles = new double[steps];
		for(int step=0; step<steps;step++) {
			angles[step] = firstAngle+angleIncr*step;
		}
		return angles;
	}
	protected static double[] getSuperpositionDistances(Atom[] ca, RotationAxis axis, double[] angles) throws StructureException {
		int steps = angles.length;

		double[] distances = new double[steps];
		if(steps < 1) return distances;

		Atom[] ca2 = StructureTools.cloneAtomArray(ca);

		// step 0
		if(angles[0] > 0) {
			axis.rotate(ca2, angles[0]);
		}
		distances[0] = superpositionDistance(ca, ca2);
		for (int step=1; step<steps;step++) {
			axis.rotate(ca2, angles[step]-angles[step-1]);
			distances[step] = superpositionDistance(ca, ca2);
		}

		return distances;
	}
	/**
	 * Fit a method-dependent function f(theta,order) to the superposition
	 * distance as atoms ca are rotated around an axis.
	 * 
	 * Returns a matrix with f(theta,order) for each angle and order in the input.
	 * Order 0 is treated as an intercept term, so that f(theta,0)=1 for all methods.
	 * 
	 * @param angles
	 * @param orders
	 * @return A matrix with angles.length rows and order.length columns.
	 */
	private Matrix computeFeatureMatrix(double[] angles, int[] orders ) {
		int steps = angles.length;
		// holds the f(angle,order) terms
		double[][] features = new double[steps][orders.length];

		for (int step=0; step<steps;step++) {
			for( int orderNum = 0; orderNum<orders.length;orderNum++) {
				int order = orders[orderNum];
				double angle = angles[step];
				double x;
				if( order == 0) {
					//initialize intercept column
					x = 1.;
				} else {
					//order-dependent column
					switch(method) {
					case HARMONICS:
					case HARMONICS_FLOATING:
					case SINGLE_HARMONIC_AMP:
					case SINGLE_HARMONIC_SSE:
						x = sin( (order)*angle/2);
						x = x*x;
						break;
					case SINGLE_CUSP_AMP:
					case SINGLE_CUSP_SSE:
						x = sqrt( 1 - cos(order*angle) );
						break;
					case SINGLE_CUSP_FIXED_AMP:
					case SINGLE_CUSP_FIXED_SSE:
						double triangleX = abs( abs(order*angle/2/PI-.5)%1 - .5);
						x = sqrt( 2 - 2*cos(2*PI/order*triangleX) );
						break;
					default:
						throw new UnsupportedOperationException("Unimplemented method "+method);
					}

				}
				features[step][orderNum] = x;
			}
		}
		return new Matrix(features);
	}

	/**
	 * Fit a method-dependent function f(theta,order) to the superposition
	 * distance as atoms ca are rotated around an axis.
	 * 
	 * Returns the root mean squared error (square root of the SSE).
	 * @param ca Atoms to rotate
	 * @param axis Axis about which to rotate ca
	 * @param orders A list of orders to include in the fit, with 0 indicating an intercept
	 * @return root sum squared error
	 * @throws StructureException For errors during rotation
	 */
	private double getSSEForFit(Atom[] ca, RotationAxis axis, int[] orders) throws StructureException {
		double[] angles = getAngles();
		double[] distances = getSuperpositionDistances(ca,axis, angles);
		return getSSEForFit(angles,distances,orders);
	}
	protected double getSSEForFit(double[] angles,double[] distances, int[] orders) throws StructureException {

		int steps = angles.length;

		Matrix y = new Matrix(distances,steps);
		Matrix M = computeFeatureMatrix(angles, orders);

		Matrix A = M.transpose().times(M).times(2);
		Matrix b = M.transpose().times(y).times(2);
		//Matrix c = y.transpose().times(y);

		// f(x) = x'Ax/2-bx+c
		// f'(x) = Ax-b = 0
		// Ax = b

		Matrix weights = A.solve(b);


		// Calculate sse
		Matrix predictions = M.times(weights);
		predictions.minusEquals(y);//errors
		Matrix sse = predictions.transpose().times(predictions);

		return sqrt(sse.get(0, 0)/steps);

	}

	/**
	 * Fit a method-dependent function f(theta,order) to the superposition
	 * distance as atoms ca are rotated around an axis.
	 * 
	 * Returns the root mean squared error (square root of the SSE).
	 * @param ca Atoms to rotate
	 * @param axis Axis about which to rotate ca
	 * @param orders A list of orders to include in the fit, with 0 indicating an intercept
	 * @return root sum squared error
	 * @throws StructureException For errors during rotation
	 */
	private double[] getWeightsForFit(Atom[] ca, RotationAxis axis, int[] orders) throws StructureException {
		double[] angles = getAngles();
		double[] distances = getSuperpositionDistances(ca,axis, angles);
		return getWeightsForFit(angles, distances, orders);
	}
	protected double[] getWeightsForFit(double[] angles, double[] distances, int[] orders) throws StructureException {

		int steps = angles.length;

		Matrix y = new Matrix(distances,steps);
		Matrix M = computeFeatureMatrix(angles, orders);

		Matrix A = M.transpose().times(M).times(2);
		Matrix b = M.transpose().times(y).times(2);
		//Matrix c = y.transpose().times(y);

		// f(x) = x'Ax/2-bx+c
		// f'(x) = Ax-b = 0
		// Ax = b

		Matrix weights = A.solve(b);

		return weights.getColumnPackedCopy();
	}


	/**
	 * For each order from 1 to maxOrder, calculate the root SSE from fitting
	 * a single-order function (with intercept).
	 * @param ca
	 * @param axis
	 * @return An array of length maxOrder containing the root mean squared error
	 * @throws StructureException For errors applying the rotation
	 */
	public double[] trySingleOrdersBySSE(Atom[] ca, RotationAxis axis) throws StructureException {
		double[] sses = new double[maxOrder];

		for( int order=1;order <= maxOrder; order++) {
			// Calculate RSSE (could save some arithmetic, but this is easier to compare)
			sses[order-1] = getSSEForFit(ca, axis, new int[] {0,order});
		}

		return sses;
	}
	/**
	 * For each order from 1 to maxOrder, calculate the amplitude from fitting
	 * a single-order function (with intercept).
	 * @param ca
	 * @param axis
	 * @return An array of length maxOrder containing the amplitude of the function
	 * @throws StructureException For errors applying the rotation
	 */
	public double[] trySingleOrdersByAmp(Atom[] ca, RotationAxis axis) throws StructureException {
		double[] amps = new double[maxOrder];

		for( int order=1;order <= maxOrder; order++) {
			// Calculate RSSE (could save some arithmetic, but this is easier to compare)
			double[] weights = getWeightsForFit(ca, axis, new int[] {0,order});
			amps[order-1] = weights[1];
		}

		return amps;
	}
	/**
	 * For each order from 1 to maxOrder, calculate the amplitude from fitting
	 * a single-order function (with intercept).
	 * @param ca
	 * @param axis
	 * @param orders array of orders to compute
	 * @return An array of length maxOrder containing the amplitude of the function
	 * @throws StructureException For errors applying the rotation
	 */
	public double[] trySingleOrdersByAmp(Atom[] ca, RotationAxis axis,int[] orders) throws StructureException {
		double[] amps = new double[orders.length];

		for( int i = 0;i<orders.length;i++) {
			// Calculate RSSE (could save some arithmetic, but this is easier to compare)
			double[] weights = getWeightsForFit(ca, axis, new int[] {0,orders[i]});
			amps[i] = weights[1];
		}

		return amps;
	}
	
	/**
	 * Fit a linear sum of f(theta,order) for all orders from 1 to maxOrder.
	 * @param ca
	 * @param axis
	 * @param intercept Indicates whether the intercept should be included (true)
	 *  or forced to 0 (false)
	 * @return An array of length maxOrder+1 containing the intercept followed
	 *  by the amplitudes for each order. If (!intercept), the first element will
	 *  always be 0.
	 * @throws StructureException 
	 */
	public double[] tryAllOrders(Atom[] ca, RotationAxis axis,boolean intercept) throws StructureException {
		if(intercept) {
			int[] orders = new int[maxOrder+1];
			for(int i=0;i<=maxOrder;i++) {
				orders[i] = i;
			}
			return getWeightsForFit(ca, axis, orders);
		} else {
			int[] orders = new int[maxOrder];
			for(int i=0;i<maxOrder;i++) {
				orders[i] = i+1;
			}
			double[] amps = getWeightsForFit(ca, axis, orders);
			// Prepend 0 for intercept
			double[] ampIntercept = new double[maxOrder+1];
			ampIntercept[0] = 0.;
			for(int i=0;i<maxOrder;i++) {
				ampIntercept[i+1] = amps[i];
			}
			return ampIntercept;
		}
	}
	
	private static int maxIndex(double[] arr, int start) {
		int maxIndex = start;
		for(int i=start+1;i<arr.length;i++) {
			if(arr[i] > arr[maxIndex])
				maxIndex = i;

		}
		return maxIndex;
	}
	private static int minIndex(double[] arr, int start) {
		int minIndex = start;
		for(int i=start+1;i<arr.length;i++) {
			if(arr[i] < arr[minIndex])
				minIndex = i;
		}
		return minIndex;
	}
}
