package org.biojava.nbio.structure.align.symm.order;

import static java.lang.Math.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Hybrid method using AngleOrderDetectorPlus for the initial detection,
 * followed by checking possible higher orders via superposition
 * @author Spencer Bliven
 *
 */
public class HybridOrderDetector extends RotationOrderDetector implements OrderDetector {
	private static Logger logger = LoggerFactory.getLogger(HybridOrderDetector.class);
	
	// list of small primes
	private static final int[] primes = new int[] {2,3,5,7,11,13,17,19,23,29,31};

	//private int maxOrder;
	private double error;
	private boolean normalizeError;
	
	private double scoreThreshold;
	
//	public HybridOrderDetector() {
//		this(8,PI,true,.85);
//	}
	public HybridOrderDetector(int maxOrder,double error,boolean normalized, double threshold) {
		super(maxOrder,RotationOrderMethod.SINGLE_CUSP_FIXED_SSE);
		//this.maxOrder = maxOrder;
		this.error = error;
		this.normalizeError = normalized;
		this.scoreThreshold = threshold;

		if(maxOrder > primes[primes.length-1]) {
			throw new IllegalArgumentException("Orders greater than "+primes[primes.length-1]+ " are not supported.");
		}
	}
	
	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca)
			throws OrderDetectionFailedException {

		try {
			RotationAxis axis = new RotationAxis(afpChain);
			double[] angles = getAngles();
			double[] distances = getSuperpositionDistances(ca,axis, angles);
			ScoreCache scores = new ScoreCache(angles,distances);

			List<Integer> compatible = compatibleOrders(afpChain, ca);

			logger.debug("Compatible orders: {}",compatible);
			
			// C1
			if(compatible.isEmpty())
				return 1;
			
			// Calculate maximum score for compatible orders
			Iterator<Integer> compatibleIt = compatible.iterator();
			int firstOrder = compatibleIt.next();
			double bestScore = scores.get(firstOrder);
			while(compatibleIt.hasNext()) {
				int order = compatibleIt.next();
				double score = scores.get(order);
				if(score > bestScore)
					bestScore = score;
			}

			int maxPrime=0;
			while(maxPrime+1 < primes.length && primes[maxPrime+1] < getMaxOrder())
				maxPrime++;

			// Find optimal order compatible with this score
			int bestOrder = -1;
			for(int order : compatible) {
				int opt = optimizeOrder(order, bestScore, 0, maxPrime, scores);
				if(opt>bestOrder)
					bestOrder = opt;
			}
			
			// Print all scores for debugging
			StringBuilder allScores = new StringBuilder();
			for(int order=1;order<=getMaxOrder();order++) {
				allScores.append(String.format("%.3f\t", scores.get(order)));
			}
			logger.debug("Scores: {}",allScores);

			return bestOrder;

		} catch (StructureException e) {
			throw new OrderDetectionFailedException(e);
		}
	}
	/**error
	 * Returns the highest order such that the score is greater than 40% of the
	 * score for startOrder.
	 * 
	 * Preconditions:
	 * -  
	 * @param startOrder Initial order
	 * @param minFactor IndmaxOrderex of the minimum element of primes to search
	 * @param maxFactor Index of the maximum element of primes to search
	 * @param scoreCache Cache of scores for each order. Should contain at least startOrder
	 * @param angles
	 * @param distances
	 * @return
	 * @throws StructureException
	 */
	private int optimizeOrder(int startOrder, double startScore,
			int minFactor, int maxFactor, final ScoreCache scoreCache) throws StructureException
	{
		logger.trace("Optimizing orders from {} with prime factors {}-{}",startOrder,primes[minFactor],primes[maxFactor]);
		final double threshold = startScore * scoreThreshold;
		
		if( minFactor > maxFactor) {
			return -1;
		}
		
		// Degenerate if <3 samples per period
		// Fails for order > 24 for standard 5deg samples
		if( 3*getAngleIncr() > 2*PI/startOrder) {
			return -1;
		}
		
		// Compute score for this order
		double score = scoreCache.get(startOrder);
		logger.trace("Score({}) = {} {}",startOrder,score, score < threshold ? "< "+ threshold : "OK");
		if(score < threshold ) {
			// Not optimal
			return -1;
		}
		
		// Still a possible optimum, so search further
		int bestOrder = startOrder;
		for(int p=minFactor;p<=maxFactor;p++) {
			int sub = optimizeOrder(startOrder*primes[p], startScore,
					p, maxFactor, scoreCache);
			if( sub > bestOrder) {
				bestOrder = sub;
			}
		}
		return bestOrder;
	}

	/**
	 * Small helper class to calculate scores as needed
	 * @author Spencer Bliven
	 *
	 */
	private class ScoreCache{
		private Map<Integer,Double> scoreCache;
		private double[] angles;
		private double[] distances;
		public ScoreCache(double[] angles, double[] distances) {
			scoreCache = new HashMap<Integer,Double>();
			this.angles = angles;
			this.distances = distances;
		}
		
		public Double get(Integer startOrder) throws StructureException {
			if( scoreCache.containsKey(startOrder)) {
				return scoreCache.get(startOrder);
			} else {
				double score = getWeightsForFit(angles,distances, new int[] {0,startOrder})[ 1 ];
				scoreCache.put(startOrder, Math.abs(score));
				return  Math.abs(score);
			}
		}
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

		for (int order = 1; order <= getMaxOrder(); order++) {
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

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return getClass().getSimpleName()+" [getMethod()=" + getMethod() +
				", error=" + error + ", scoreThreshold=" + scoreThreshold
				+ ", getMaxOrder()=" + getMaxOrder() + ", getAngleIncr()="
				+ getAngleIncr() + "]";
	}

}
