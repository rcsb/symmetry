package org.biojava.nbio.structure.align.symm.census3.stats.order;

import java.util.HashMap;
import java.util.Map;

import org.biojava.nbio.structure.align.symm.census3.CensusResult;
import org.biojava.nbio.structure.align.symm.census3.stats.CensusStatUtils;
import org.biojava.nbio.structure.align.symm.census3.stats.StructureClassificationGrouping;

/**
 * A helper class that can decide the order of a fold based on consensus of the domains.
 * 
 * @author dmyersturnbull
 */
public class OrderHelper {

	private Map<String, Integer> hasOrderMap = new HashMap<String, Integer>();
	private Map<String, Map<Integer,Integer>> ordersMap = new HashMap<String, Map<Integer,Integer>>();
	private Map<String, Float> tmScoreMap = new HashMap<String, Float>();
	private Map<String, Integer> counts = new HashMap<String, Integer>();

	private StructureClassificationGrouping normalizer;
	private double tmScoreCutoff;
	private ConsensusDecider consensusDecider;
	
	public OrderHelper(StructureClassificationGrouping normalizer, double tmScoreCutoff, ConsensusDecider consensusDecider) {
		super();
		this.normalizer = normalizer;
		this.tmScoreCutoff = tmScoreCutoff;
		this.consensusDecider = consensusDecider;
	}

	public void add(CensusResult result) {

		if (result.getAlignment() == null) return;

		String fold = normalizer.group(result);

		CensusStatUtils.plus(counts, fold);

		CensusStatUtils.plusF(tmScoreMap, fold, result.getScoreList().getTmScore());

		// also note it would be unfair to decrease the hasOrder score with -1
		Integer order = result.getOrder();
		if (order == null || result.getOrder() < 2) order = 1;
		CensusStatUtils.plus(hasOrderMap, fold, order > 1 ? 1 : 0);

		// do this AFTER we do hasOrderMap, but before ordersMap
//		if (!sig.isSignificant(result)) order = 1;
		
		if (ordersMap.get(fold) == null) {
			ordersMap.put(fold, new HashMap<Integer,Integer>(7));
		}
		CensusStatUtils.plus(ordersMap.get(fold), order);

	}

	/**
	 * Returns 1 if, among domains in {@code fold}:
	 * <ul>
	 * <li>The mean TM-score &lt; tmScoreCutoff</li>
	 * <li>The mean [order>1] &lt; orderScoreCutoff</li>
	 * </ul>
	 * Otherwise, returns the mode of orders among symmetric domains. Here, significance is decided by the criteria
	 * [order &gt; 1 and TM-score at least tmScoreCutoff].
	 * 
	 * @param fold
	 *            A SCOP fold, superfamily, or family given by {@code normalizer} in the parent class
	 */
	public int getConsensusOrder(String fold) {

		// if asymmetric, return 1
		if (!counts.containsKey(fold)) return 1;
		int n = counts.get(fold);
		Float tmScore = tmScoreMap.get(fold);
		if (tmScore == null) tmScore = 0.0f;
		Integer order = hasOrderMap.get(fold);
		if (order == null) order = 1;
		double meanTmScore = (double) tmScore / n;
		if (meanTmScore < tmScoreCutoff) {
			return 1;
		}

		// otherwise return the consensus order
		int maximizingOrder = consensusDecider.decide(ordersMap.get(fold));
		assert(maximizingOrder > 0);
		return maximizingOrder;
	}

}
