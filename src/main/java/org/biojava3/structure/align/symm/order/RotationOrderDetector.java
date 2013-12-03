package org.biojava3.structure.align.symm.order;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JOptionPane;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.align.symm.CeSymm;

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

	private static final double DEFAULT_ANGLE_INCR = 5*Calc.radiansPerDegree;
	private final double angleIncr = DEFAULT_ANGLE_INCR; // angular resolution

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
	public static void printSuperpositionDistance(Atom[] ca, RotationAxis axis, PrintStream out) throws StructureException {
		printSuperpositionDistance(ca, axis, DEFAULT_ANGLE_INCR, out);
	}
	public static void printSuperpositionDistance(Atom[] ca, RotationAxis axis,double angleIncr, PrintStream out) throws StructureException {
		final int steps = (int)Math.floor(Math.PI/angleIncr);

		Atom[] ca2 = StructureTools.cloneCAArray(ca);

		for (int step=0; step<steps;step++) {
			double dist = superpositionDistance(ca, ca2);
			double angle = angleIncr*step;

			// Rotate for next step
			axis.rotate(ca2, angleIncr);

			out.format("%f\t%f%n",angle,dist);
		}
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

	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		//argument parsing
		final String usage = "[OPTIONS] [structure]";
		final String header = "Determine the order for <structure>, which may " +
				"be a PDB ID, SCOP domain, or file path. If none is given, the " +
				"user will be prompted at startup.";
		Options options = new Options();
		options.addOption("h","help",false,"print help");
		options.addOption(OptionBuilder.withArgName("maxorder")
				.hasArg()
				.withLongOpt("max-order")
				.withType(Number.class)
				.withDescription("maximum order to consider [default 8]")
				.create("x") );
		StringBuilder descr = new StringBuilder("Method to use. May be specified multiple times. [");
		RotationOrderMethod[] availableMethods = RotationOrderMethod.values();
		descr.append(availableMethods[0]);
		for(int i=1;i<availableMethods.length;i++) {
			descr.append(" | ");
			descr.append(availableMethods[i]);
		}
		descr.append(']');
		options.addOption(OptionBuilder.withArgName("method")
				.hasArg()
				.withLongOpt("method")
				.withType(RotationOrderMethod.class)
				.withDescription(descr.toString())
				.create('M') );
		options.addOption("o","output",true,"tab delimited output file containing angle and distance columns");
		CommandLineParser parser = new GnuParser();
		HelpFormatter help = new HelpFormatter();

		CommandLine cli;
		try {
			cli = parser.parse(options,args,false);
			if(cli.hasOption('h')) {
				help.printHelp(usage, header, options, "");
				System.exit(1);
				return;
			}
		} catch (ParseException e) {
			System.err.println("Error: "+e.getMessage());
			help.printHelp(usage, header, options, "");
			System.exit(1);
			return;
		}

		args = cli.getArgs();


		String name;
		if(args.length == 0) {
			// default name
			name = "d1ijqa1";
			//		name = "1G6S";
			name = "1MER.A";
			//		name = "1MER";
			//		name = "1TIM.A";
			//		name = "d1h70a_";
			name = "2YMS";

			name = (String)JOptionPane.showInputDialog(
					null,
					"Structure ID (PDB, SCOP, etc):",
					"Input Structure",
					JOptionPane.PLAIN_MESSAGE,
					null,
					null,
					name);

			if( name == null) {
				//cancel
				return;
			}
		} else if(args.length == 1) {
			name = args[0];
		} else {
			help.printHelp(usage, header, options, "");
			System.exit(1);
			return;
		}

		int maxorder = 8;
		if(cli.hasOption('x') ) {
			maxorder = Integer.parseInt(cli.getOptionValue('x'));
		}

		String outfile = null;
		if(cli.hasOption('o')) {
			outfile = cli.getOptionValue('o');
		}

		List<RotationOrderDetector> methods = new ArrayList<RotationOrderDetector>();
		if(cli.hasOption('M')) {
			for(String method: cli.getOptionValues('M')) {
				RotationOrderMethod m = RotationOrderMethod.valueOf(method);
				methods.add(new RotationOrderDetector(maxorder,m) );
			}
		}

//		System.out.println("Name:" + name);
//		System.out.println("order:" + maxorder);
//		System.out.println("output:" + outfile);
//		for(RotationOrderDetector m: methods) {
//			System.out.println("method:"+m);
//		}

		// Done parsing arguments


		try {

			// Perform alignment to determine axis
			Atom[] ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
			Atom[] ca2 = StructureTools.cloneCAArray(ca1);
			CeSymm ce = new CeSymm();
			AFPChain alignment = ce.align(ca1, ca2);
			RotationAxis axis = new RotationAxis(alignment);

			// Output raw data
			if(outfile != null) {
				PrintStream out = null;
				try {
					out = new PrintStream(outfile);
					out.println("Angle\tDistance");
					printSuperpositionDistance(ca1,axis,out);
				} catch(FileNotFoundException e) {
					e.printStackTrace();
				} finally {
					if(out != null) {
						out.close();
					}
				}
			}

			// Print orders
			for( OrderDetector detector: methods) {
				// Calculate order
				int order = detector.calculateOrder(alignment, ca1);
				System.out.format("%s\t%d%n",detector,order);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (StructureException e) {
			e.printStackTrace();
		} catch (OrderDetectionFailedException e) {
			e.printStackTrace();
		}



	}
}
