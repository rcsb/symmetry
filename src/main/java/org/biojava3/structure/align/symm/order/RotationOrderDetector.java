package org.biojava3.structure.align.symm.order;

import java.io.IOException;
import java.util.Arrays;

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

	private int maxOrder;

	public RotationOrderDetector() {
		this(8);
	}

	public RotationOrderDetector(int maxOrder) {
		super();
		this.maxOrder = maxOrder;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {

		try {

			RotationAxis axis = new RotationAxis(afpChain);

			// Use C1 order if the axis is undefined
			if(!axis.isDefined()) {
				return 1;
			}
			// Calculate weights for each order
			//TODO only use aligned residues, rather than the whole ca
			double[] coefficients = fitHarmonics(ca,axis);

			// Find order with maximum weight
			double bestScore = coefficients[0];
			int argmax = 1;
			for(int i=1;i<coefficients.length;i++) {
				if(coefficients[i] > bestScore ) {
					argmax = i+1;
					bestScore = coefficients[i];
				}
			}
			return argmax;

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
	public double superpositionDistance(Atom[] ca1, Atom[] ca2) throws StructureException {

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
	 *  f(angle) = a1*sin(angle/2)^2 + a2*sin(2*angle/2)^2 + a3*sin(3*angle/2)^2 + ...
	 *  
	 * Performs a best-fit to find coefficients a1,a2,...,aMaxOrder.
	 * These give the relative contribution of each order to the overall function.
	 * @param ca Aligned residues from the protein structure
	 * @param axis The rotaton axis about which to rotate
	 * @return an array of length maxOrder giving the coefficients
	 * @throws StructureException
	 */
	public double[] fitHarmonics(Atom[] ca, RotationAxis axis) throws StructureException {
		final double angleIncr = 5*Calc.radiansPerDegree;
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
		//y=y.getMatrix(1, steps-1, 0, 0);
		//M=M.getMatrix(1, steps-1, 1, maxOrder-1);
		
		
		Matrix A = M.transpose().times(M).times(2);
		Matrix b = M.transpose().times(y).times(2);
		Matrix c = y.transpose().times(y);
		
		// f(x) = x'Ax/2-bx+c
		// f'(x) = Ax-b = 0
		// Ax = b
//		assert(A.getRowDimension() == maxOrder);
//		assert(A.getColumnDimension() == maxOrder);
//		assert(b.getRowDimension() == 1);
//		assert(b.getColumnDimension() == maxOrder);
//		assert(c.getRowDimension() == 1);
//		assert(c.getColumnDimension() == 1);
		
		A.inverse();
		
		Matrix harmonicWeights = A.solve(b);
		
//		assert(harmonicWeights.getRowDimension() == maxOrder);
//		assert(harmonicWeights.getColumnDimension() == 1);

		return harmonicWeights.transpose().getArray()[0];
	}
	
	public static void main(String[] args) {
		String name;
		name = "d1ijqa1";
//		name = "1G6S";
		name = "1MER.A";
//		name = "1MER";
//		name = "1TIM.A";
//		name = "d1h70a_";
		name = "2YMS";
		try {
			
			// Perform alignment to determine axis
			Atom[] ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
			Atom[] ca2 = StructureTools.cloneCAArray(ca1);
			CeSymm ce = new CeSymm();
			AFPChain alignment = ce.align(ca1, ca2);

			// Search for orders up to 9
			RotationOrderDetector detector = new RotationOrderDetector(9);
			
			// Calculate order
			int order = detector.calculateOrder(alignment, ca1);
			System.out.println("Order: "+order);

			// Print weights for each order, for comparison
			RotationAxis axis = new RotationAxis(alignment);
			double[] harmonics = detector.fitHarmonics(ca1,axis);
			System.out.println("Order Weights: "+Arrays.toString(harmonics));
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (StructureException e) {
			e.printStackTrace();
		} catch (OrderDetectionFailedException e) {
			e.printStackTrace();
		}
		
		
		
	}
}
