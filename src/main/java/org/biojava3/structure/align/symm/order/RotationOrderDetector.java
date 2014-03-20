package org.biojava3.structure.align.symm.order;

import java.util.Arrays;

import org.apache.commons.math3.util.Pair;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.jama.Matrix;

/**
 * Detects order by analyzing the goodness of fit as the protein is rotated
 * around the axis of symmetry.
 * @author Spencer Bliven
 */
public class RotationOrderDetector implements OrderDetector {
	public static enum RotationOrderMethod {
		HARMONICS,
		HARMONICS_FLOATING,
		SINGLE_HARMONIC_AMP,
		SINGLE_HARMONIC_SSE,
		SINGLE_CUSP_AMP,
		SINGLE_CUSP_SSE
	}

	private int maxOrder;

	public static final double DEFAULT_ANGLE_INCR = 5*Calc.radiansPerDegree;
	private double angleIncr = DEFAULT_ANGLE_INCR; // angular resolution

	private RotationOrderMethod method;
	public RotationOrderDetector() {
		this(8);
	}

	public RotationOrderDetector(int maxOrder) {
		this(maxOrder,RotationOrderMethod.SINGLE_HARMONIC_SSE);
	}

	public RotationOrderDetector(int maxOrder, RotationOrderMethod method) {
		this.maxOrder = maxOrder;
		this.method= method;
	}

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
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {
		//TODO only use aligned residues, rather than the whole ca
		switch(method) {
		case HARMONICS:
			return calculateOrderHarmonics(afpChain,ca);
		case HARMONICS_FLOATING:
			return calculateOrderHarmonicsFloating(afpChain,ca);
		case SINGLE_HARMONIC_AMP:
			return calculateOrderSingleHarmonicsByAmp(afpChain,ca);
		case SINGLE_HARMONIC_SSE:
			return calculateOrderSingleHarmonicsBySSE(afpChain,ca);
		case SINGLE_CUSP_AMP:
			return calculateOrderSingleCuspByAmp(afpChain,ca);
		case SINGLE_CUSP_SSE:
		default:
			return calculateOrderSingleCuspBySSE(afpChain,ca);
		}
	}

	/**
	 * Returns an array of {@link #superpositionDistance(Atom[], Atom[]) superposition distances} of rotations of {@code ca}.
	 * The {@code n}th element in the array corresponds to a rotation by {@code degreesIncrement * n} degrees.
	 */
	public static Pair<double[],double[]> sampleRotations(Atom[] ca, RotationAxis axis, double degreesIncrement) throws StructureException {
		
		final double angleIncr = degreesIncrement*Calc.radiansPerDegree;
		final int steps = (int)Math.floor(2*Math.PI/angleIncr);

		double[] angles = new double[steps];
		double[] distances = new double[steps];

		Atom[] ca2 = StructureTools.cloneCAArray(ca);
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
	
	int calculateOrderHarmonics(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);

			// Use C1 order if the axis is undefined
			if(!axis.isDefined()) {
				return 1;
			}
			// Calculate weights for each order
			double[] coefficients = fitHarmonics(ca,axis);

			// Find order with maximum weight
			// ignore initial intercept term
			double bestScore = coefficients[1];
			int maxorder = 1;
			for(int order=2;order<coefficients.length;order++) {
				if(coefficients[order] > bestScore ) {
					maxorder = order;
					bestScore = coefficients[order];
				}
			}
			return maxorder;

		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
	}
	int calculateOrderHarmonicsFloating(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);

			// Use C1 order if the axis is undefined
			if(!axis.isDefined()) {
				return 1;
			}
			// Calculate weights for each order
			double[] coefficients = fitHarmonicsFloating(ca,axis);

			// Find order with maximum weight
			// ignore initial intercept term
			double bestScore = coefficients[1];
			int maxorder = 1;
			for(int order=2;order<coefficients.length;order++) {
				if(coefficients[order] > bestScore ) {
					maxorder = order;
					bestScore = coefficients[order];
				}
			}
			return maxorder;

		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
	}
	int calculateOrderSingleHarmonicsByAmp(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);

			// Use C1 order if the axis is undefined
			if(!axis.isDefined()) {
				return 1;
			}
			// Calculate weights for each order
			double[] coefficients = trySingleHarmonicsFloatingByAmp(ca,axis);

			// Find order with maximum weight
			// ignore initial intercept term
			double bestScore = coefficients[0];
			int maxorder = 1;
			for(int order=2;order<coefficients.length;order++) {
				if(coefficients[order-1] > bestScore ) {
					maxorder = order;
					bestScore = coefficients[order-1];
				}
			}
			return maxorder;

		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
	}
	int calculateOrderSingleHarmonicsBySSE(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);

			// Use C1 order if the axis is undefined
			if(!axis.isDefined()) {
				return 1;
			}
			// Calculate weights for each order
			double[] coefficients = trySingleHarmonicsFloatingBySSE(ca,axis);

			// Find order with maximum weight
			// ignore initial intercept term
			double bestScore = coefficients[0];
			int maxorder = 1;
			for(int order=2;order<coefficients.length;order++) {
				if(coefficients[order-1] < bestScore ) {
					maxorder = order;
					bestScore = coefficients[order-1];
				}
			}
			return maxorder;

		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
	}
	int calculateOrderSingleCuspByAmp(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);

			// Use C1 order if the axis is undefined
			if(!axis.isDefined()) {
				return 1;
			}
			// Calculate weights for each order
			double[] coefficients = trySingleCuspByAmp(ca,axis);

			// Find order with maximum weight
			// ignore initial intercept term
			double bestScore = coefficients[0];
			int maxorder = 1;
			for(int order=2;order<coefficients.length;order++) {
				if(coefficients[order-1] > bestScore ) {
					maxorder = order;
					bestScore = coefficients[order-1];
				}
			}
			return maxorder;

		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
	}
	int calculateOrderSingleCuspBySSE(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);

			// Use C1 order if the axis is undefined
			if(!axis.isDefined()) {
				return 1;
			}
			// Calculate weights for each order
			double[] coefficients = trySingleCuspBySSE(ca,axis);

			// Find order with maximum weight
			// ignore initial intercept term
			double bestScore = coefficients[0];
			int maxorder = 1;
			for(int order=2;order<coefficients.length;order++) {
				if(coefficients[order-1] < bestScore ) {
					maxorder = order;
					bestScore = coefficients[order-1];
				}
			}
			return maxorder;

		} catch (Exception e) {
			throw new OrderDetectionFailedException(e);
		}
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
			total += Math.sqrt(bestDist1[i]);
		}
		for(int j=0;j<ca2.length;j++) {
			total += Math.sqrt(bestDist2[j]);
		}

		double dist = total/(ca1.length+ca2.length);
		return dist;
	}

	/**
	 * Models the relationship between order and superpositionDistance as a
	 * combination of sine squared functions:
	 * 
	 *  f(angle) = 0 + a1*sin(angle/2)^2 + a2*sin(2*angle/2)^2 + a3*sin(3*angle/2)^2 + ...
	 *  
	 * Performs a best-fit to find coefficients a1,a2,...,aMaxOrder.
	 * These give the relative contribution of each order to the overall function.
	 * @param ca Aligned residues from the protein structure
	 * @param axis The rotaton axis about which to rotate
	 * @return an array of length maxOrder+1 giving the intercept (0) and coefficients a1...a_maxOrder
	 * @throws StructureException
	 */
	double[] fitHarmonics(Atom[] ca, RotationAxis axis) throws StructureException {
		final int steps = (int)Math.floor(Math.PI/angleIncr);

		// Fit (angle,distance) points to the following equation
		// f(angle) = a1*sin^2(1*angle/2) + a2*sin^2(2*angle/2) +...+ a*sin2(maxOrder*angle/2)
		// goal is to find a_1...a_maxOrder

		double[][] distances = new double[steps][1];
		// holds the sin2(i*angle/2) terms
		double[][] harmonics = new double[steps][maxOrder];

		Atom[] ca2 = StructureTools.cloneCAArray(ca);

		for (int step=0; step<steps;step++) {
			double dist = superpositionDistance(ca, ca2);
			distances[step][0] = dist;

			for(int order=0;order<maxOrder;order++) {
				double angle = angleIncr*step;
				double x = Math.sin( (order+1)*angle/2);
				harmonics[step][order] = x*x;
			}
			// Rotate for next step
			axis.rotate(ca2, angleIncr);
		}

		Matrix y = new Matrix(distances);
		Matrix M = new Matrix(harmonics);

		Matrix A = M.transpose().times(M).times(2);
		Matrix b = M.transpose().times(y).times(2);
		//Matrix c = y.transpose().times(y);

		// f(x) = x'Ax/2-bx+c
		// f'(x) = Ax-b = 0
		// Ax = b

		Matrix harmonicWeights = A.solve(b);

		double[] coefs = new double[maxOrder+1];
		coefs[0] = 0.;
		for(int i=0;i<maxOrder;i++) {
			coefs[i+1] = harmonicWeights.get( i, 0);
		}
		return coefs;
	}

	/**
	 * Like fitHarmonics, but adds an intercept term and ignores points too close
	 * to angle 0, which artificially is forced to zero.
	 * 
	 *  f(angle) = a0 + a1*sin(angle/2)^2 + a2*sin(2*angle/2)^2 + a3*sin(3*angle/2)^2 + ...
	 *  
	 * Performs a best-fit to find coefficients a1,a2,...,aMaxOrder.
	 * These give the relative contribution of each order to the overall function.
	 * @param ca Aligned residues from the protein structure
	 * @param axis The rotaton axis about which to rotate
	 * @return an array of length maxOrder+1 giving the coefficients a0...a_maxOrder
	 *  (intercept is element 0)
	 * @throws StructureException
	 */
	double[] fitHarmonicsFloating(Atom[] ca, RotationAxis axis) throws StructureException {
		// Range of angles to use for training
		final double minAngle = Math.floor(Math.PI/maxOrder/angleIncr)*angleIncr; // first valid peak
		final double maxAngle = Math.PI;
		// Number of angle steps
		final int steps = (int)Math.floor((maxAngle-minAngle)/angleIncr);

		// Fit (angle,distance) points to the following equation
		// f(angle) = a0 + a1*sin^2(1*angle/2) + a2*sin^2(2*angle/2) +...+ a*sin2(maxOrder*angle/2)
		// goal is to find a_1...a_maxOrder

		double[][] distances = new double[steps][1];//preserve matrix dimensions
		// holds the sin2(i*angle/2) terms
		double[][] harmonics = new double[steps][maxOrder+1];

		Atom[] ca2 = StructureTools.cloneCAArray(ca);

		if(minAngle != 0.) {
			axis.rotate(ca2, minAngle);
		}
		for (int step=0; step<steps;step++) {
			double dist = superpositionDistance(ca, ca2);
			distances[step][0] = dist;

			//initialize intercept column
			harmonics[step][0] = 1.;

			//later columns
			for(int order=1;order<=maxOrder;order++) {
				double angle = minAngle+angleIncr*step;
				double x = Math.sin( order*angle/2);
				harmonics[step][order] = x*x;
			}
			// Rotate for next step
			axis.rotate(ca2, angleIncr);
		}

		Matrix y = new Matrix(distances);
		Matrix M = new Matrix(harmonics);

		Matrix A = M.transpose().times(M).times(2);
		Matrix b = M.transpose().times(y).times(2);
		//Matrix c = y.transpose().times(y);

		// f(x) = x'Ax/2-bx+c
		// f'(x) = Ax-b = 0
		// Ax = b

		Matrix harmonicWeights = A.solve(b);

		return harmonicWeights.transpose().getArray()[0];
	}

	/**
	 * Like fitHarmonics, but adds an intercept term and ignores points too close
	 * to angle 0, which artificially is forced to zero.
	 * 
	 *  f(angle) = a0 + a1*sin(angle/2)^2 + a2*sin(2*angle/2)^2 + a3*sin(3*angle/2)^2 + ...
	 *  
	 * Performs a best-fit to find coefficients a1,a2,...,aMaxOrder.
	 * These give the relative contribution of each order to the overall function.
	 * @param ca Aligned residues from the protein structure
	 * @param axis The rotaton axis about which to rotate
	 * @return an array of length maxOrder+1 giving the coefficients a0...a_maxOrder
	 *  (intercept is element 0)
	 * @throws StructureException
	 */
	double[] trySingleHarmonicsFloatingByAmp(Atom[] ca, RotationAxis axis) throws StructureException {
		// Range of angles to use for training
		final double minAngle = Math.floor(Math.PI/maxOrder/angleIncr)*angleIncr; // first valid peak
		final double maxAngle = Math.PI;
		// Number of angle steps
		final int steps = (int)Math.floor((maxAngle-minAngle)/angleIncr);

		// Fit (angle,distance) points to the following equation
		// f(angle) = a0 + a1*sin^2(1*angle/2) + a2*sin^2(2*angle/2) +...+ a*sin2(maxOrder*angle/2)
		// goal is to find a_1...a_maxOrder

		double[][] distances = new double[steps][1];//preserve matrix dimensions
		Atom[] ca2 = StructureTools.cloneCAArray(ca);
		if(minAngle != 0.) {
			axis.rotate(ca2, minAngle);
		}
		for (int step=0; step<steps;step++) {
			double dist = superpositionDistance(ca, ca2);
			distances[step][0] = dist;
			// Rotate for next step
			axis.rotate(ca2, angleIncr);
		}

		double[] amplitudes = new double[maxOrder];

		for( int order=1;order <= maxOrder; order++) {

			// holds the sin2(i*angle/2) terms
			double[][] harmonics = new double[steps][2];

			for (int step=0; step<steps;step++) {

				//initialize intercept column
				harmonics[step][0] = 1.;

				//order-dependent column
				double angle = minAngle+angleIncr*step;
				double x = Math.sin( order*angle/2);
				harmonics[step][1] = x*x;
			}

			Matrix y = new Matrix(distances);
			Matrix M = new Matrix(harmonics);

			Matrix A = M.transpose().times(M).times(2);
			Matrix b = M.transpose().times(y).times(2);
			//Matrix c = y.transpose().times(y);

			// f(x) = x'Ax/2-bx+c
			// f'(x) = Ax-b = 0
			// Ax = b

			Matrix harmonicWeights = A.solve(b);

			amplitudes[order-1] = harmonicWeights.get(1, 0);
		}

		return amplitudes;
		//return sses;

		//		int bestOrder = 1;
		//		double minScore = amplitudes[0];
		//		for(int order=2;order<maxOrder;order++) {
		//			if(amplitudes[order-1] < minScore) {
		//				minScore = amplitudes[order-1];
		//				bestOrder = order;
		//			}
		//		}
		//		return bestOrder;
	}

	double[] trySingleHarmonicsFloatingBySSE(Atom[] ca, RotationAxis axis) throws StructureException {
		// Range of angles to use for training
		final double minAngle = Math.floor(Math.PI/maxOrder/angleIncr)*angleIncr; // first valid peak
		final double maxAngle = Math.PI;
		// Number of angle steps
		final int steps = (int)Math.floor((maxAngle-minAngle)/angleIncr);

		// Fit (angle,distance) points to the following equation
		// f(angle) = a0 + a1*sin^2(1*angle/2) + a2*sin^2(2*angle/2) +...+ a*sin2(maxOrder*angle/2)
		// goal is to find a_1...a_maxOrder

		double[][] distances = new double[steps][1];//preserve matrix dimensions
		Atom[] ca2 = StructureTools.cloneCAArray(ca);
		if(minAngle != 0.) {
			axis.rotate(ca2, minAngle);
		}
		for (int step=0; step<steps;step++) {
			double dist = superpositionDistance(ca, ca2);
			distances[step][0] = dist;
			// Rotate for next step
			axis.rotate(ca2, angleIncr);
		}

		double[] sses = new double[maxOrder];

		for( int order=1;order <= maxOrder; order++) {

			// holds the sin2(i*angle/2) terms
			double[][] harmonics = new double[steps][2];

			for (int step=0; step<steps;step++) {

				//initialize intercept column
				harmonics[step][0] = 1.;

				//order-dependent column
				double angle = minAngle+angleIncr*step;
				double x = Math.sin( order*angle/2);
				harmonics[step][1] = x*x;
			}

			Matrix y = new Matrix(distances);
			Matrix M = new Matrix(harmonics);

			Matrix A = M.transpose().times(M).times(2);
			Matrix b = M.transpose().times(y).times(2);
			//Matrix c = y.transpose().times(y);

			// f(x) = x'Ax/2-bx+c
			// f'(x) = Ax-b = 0
			// Ax = b

			Matrix harmonicWeights = A.solve(b);

			Matrix predictions = M.times(harmonicWeights);
			predictions.minusEquals(y);//errors
			Matrix sse = predictions.transpose().times(predictions);

			// Calculate RSSE (could save some arithmetic, but this is easier to compare)
			sses[order-1] = Math.sqrt(sse.get(0, 0)/steps);
		}

		return sses;
	}
	double[] trySingleCuspByAmp(Atom[] ca, RotationAxis axis) throws StructureException {
		// Range of angles to use for training
		final double minAngle = Math.floor(Math.PI/maxOrder/angleIncr)*angleIncr; // first valid peak
		final double maxAngle = Math.PI;
		// Number of angle steps
		final int steps = (int)Math.floor((maxAngle-minAngle)/angleIncr);

		// Fit (angle,distance) points to the following equation
		// f(angle) = a0 + a1*sin^2(1*angle/2) + a2*sin^2(2*angle/2) +...+ a*sin2(maxOrder*angle/2)
		// goal is to find a_1...a_maxOrder

		double[][] distances = new double[steps][1];//preserve matrix dimensions
		Atom[] ca2 = StructureTools.cloneCAArray(ca);
		if(minAngle != 0.) {
			axis.rotate(ca2, minAngle);
		}
		for (int step=0; step<steps;step++) {
			double dist = superpositionDistance(ca, ca2);
			distances[step][0] = dist;
			// Rotate for next step
			axis.rotate(ca2, angleIncr);
		}

		double[] amplitudes = new double[maxOrder];

		for( int order=1;order <= maxOrder; order++) {

			// holds the sin2(i*angle/2) terms
			double[][] harmonics = new double[steps][2];

			for (int step=0; step<steps;step++) {

				//initialize intercept column
				harmonics[step][0] = 1.;

				//order-dependent column
				double angle = minAngle+angleIncr*step;
				double x = Math.sqrt( 1 - Math.cos(order*angle) );
				harmonics[step][1] = x;
			}

			Matrix y = new Matrix(distances);
			Matrix M = new Matrix(harmonics);

			Matrix A = M.transpose().times(M).times(2);
			Matrix b = M.transpose().times(y).times(2);
			//Matrix c = y.transpose().times(y);

			// f(x) = x'Ax/2-bx+c
			// f'(x) = Ax-b = 0
			// Ax = b

			Matrix harmonicWeights = A.solve(b);

			amplitudes[order-1] = harmonicWeights.get(1, 0);
		}

		return amplitudes;
	}
	double[] trySingleCuspBySSE(Atom[] ca, RotationAxis axis) throws StructureException {
		// Range of angles to use for training
		final double minAngle = Math.floor(Math.PI/maxOrder/angleIncr)*angleIncr; // first valid peak
		final double maxAngle = Math.PI;
		// Number of angle steps
		final int steps = (int)Math.floor((maxAngle-minAngle)/angleIncr);

		// Fit (angle,distance) points to the following equation
		// f(angle) = a0 + a1*sin^2(1*angle/2) + a2*sin^2(2*angle/2) +...+ a*sin2(maxOrder*angle/2)
		// goal is to find a_1...a_maxOrder

		double[][] distances = new double[steps][1];//preserve matrix dimensions
		Atom[] ca2 = StructureTools.cloneCAArray(ca);
		if(minAngle != 0.) {
			axis.rotate(ca2, minAngle);
		}
		for (int step=0; step<steps;step++) {
			double dist = superpositionDistance(ca, ca2);
			distances[step][0] = dist;
			// Rotate for next step
			axis.rotate(ca2, angleIncr);
		}

		double[] sses = new double[maxOrder];

		for( int order=1;order <= maxOrder; order++) {

			// holds the sin2(i*angle/2) terms
			double[][] harmonics = new double[steps][2];

			for (int step=0; step<steps;step++) {

				//initialize intercept column
				harmonics[step][0] = 1.;

				//order-dependent column
				double angle = minAngle+angleIncr*step;
				double x = Math.sqrt( 1 - Math.cos(order*angle) );
				harmonics[step][1] = x*x;
			}

			Matrix y = new Matrix(distances);
			Matrix M = new Matrix(harmonics);

			Matrix A = M.transpose().times(M).times(2);
			Matrix b = M.transpose().times(y).times(2);
			//Matrix c = y.transpose().times(y);

			// f(x) = x'Ax/2-bx+c
			// f'(x) = Ax-b = 0
			// Ax = b

			Matrix harmonicWeights = A.solve(b);

			Matrix predictions = M.times(harmonicWeights);
			predictions.minusEquals(y);//errors
			Matrix sse = predictions.transpose().times(predictions);

			// Calculate RSSE (could save some arithmetic, but this is easier to compare)
			sses[order-1] = Math.sqrt(sse.get(0, 0)/steps);
		}

		return sses;
	}
}
