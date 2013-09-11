package org.biojava3.structure.align.symm.census2.stats;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;

/**
 * Tabulate symmetry by order, where symmetry order is determined by consensus over domains in a SCOP category.
 * 
 * @author dmyerstu
 */
public class FoldOrder {

	public interface ConsensusDecider {
		int decide(Map<Integer,Integer> orders);
	}

	public static final ConsensusDecider MODE_DECIDER = new ConsensusDecider() {
		@Override
		public int decide(Map<Integer,Integer> countsByOrders) {
			int maxCount = 0;
			int maximizingOrder = 0;
			for (int i = 2; i <= 8; i++) {
				Integer count = countsByOrders.get(i);
				if (count == null) count = 0;
				if (count > maxCount) {
					maxCount = count;
					maximizingOrder = i;
				}
			}
			return maximizingOrder;
		}
	};

	/**
	 * A {@link CensusDecider} that treats any <em>divisor</em> of an order as potentially being the order itself.
	 * TODO: Decide increases with a little more basis.
	 */
	public static final ConsensusDecider SMART_DECIDER = new ConsensusDecider() {
		@Override
		public int decide(Map<Integer,Integer> countsByOrders) {
			int[] newCounts = new int[] {0, 0, 0, 0, 0, 0, 0};
			for (int i = 2; i <= 8; i++) {
				for (int j = 2; j <= 8; j++) {
					/*
					 * Ex: If i=6 and j=3, say that each j is the square root of an i
					 */
					Integer c = countsByOrders.get(j);
					if (c == null) c = 0;
					double increase = Math.pow((double) c, (double) j / (double) i);
//					for (int p = 1; p < i / j; p++) { // make sure we don't do this loop if i = j
//						increase = Math.log(increase);
//					}
					if (i % j == 0) newCounts[i-2] += increase;
				}
			}
			int maxCount = 0;
			int maximizingOrder = 0;
			for (int i = 2; i <= 8; i++) {
				int count = newCounts[i-2];
				if (count > maxCount) {
					maxCount = count;
					maximizingOrder = i;
				}
			}
			return maximizingOrder;
		}
	};
	
	private ConsensusDecider consensusDecider = SMART_DECIDER;
	
	/**
	 * A helper class that can decide the order of a fold based on consensus of the domains.
	 * 
	 * @author dmyerstu
	 */
	public class OrderDecider {

		private Map<String, Integer> hasOrderMap = new HashMap<String, Integer>();
		private Map<String, Map<Integer,Integer>> ordersMap = new HashMap<String, Map<Integer,Integer>>();
		private Map<String, Float> tmScoreMap = new HashMap<String, Float>();
		private Map<String, Integer> counts = new HashMap<String, Integer>();

		public void add(Result result) {

			if (result.getAlignment() == null) return;

			String fold = normalizer.group(result);

			StatUtils.plus(counts, fold);

			StatUtils.plusF(tmScoreMap, fold, result.getAlignment().getTmScore());

			// also note it would be unfair to decrease the hasOrder score with -1
			Integer order = result.getOrder();
			if (order == null || result.getOrder() < 2) order = 1;
			StatUtils.plus(hasOrderMap, fold, order > 1 ? 1 : 0);

			if (ordersMap.get(fold) == null) {
				ordersMap.put(fold, new HashMap<Integer,Integer>(7));
			}
			StatUtils.plus(ordersMap.get(fold), order);

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
			int n = counts.get(fold);
			Float tmScore = tmScoreMap.get(fold);
			if (tmScore == null) tmScore = 0.0f;
			Integer order = hasOrderMap.get(fold);
			if (order == null) order = 1;
			double meanTmScore = (double) tmScore / n;
			double meanHasOrder = (double) order / n;
			if (meanTmScore < tmScoreCutoff || meanHasOrder < orderCutoff) {
				return 1;
			}

			// otherwise return the consensus order
			int maximizingOrder = consensusDecider.decide(ordersMap.get(fold));
			assert(maximizingOrder > 0);
			return maximizingOrder;
		}

	}

	private static final Logger logger = LogManager.getLogger(SymmetryOrder.class.getName());

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2) {
			System.err.println("Usage: " + SymmetryOrder.class.getSimpleName() + " census-file.xml [grouping]");
			return;
		}
		FoldOrder orders = new FoldOrder();
		if (args.length > 1) {
			orders.setNormalizer(Grouping.byName(args[1]));
		}
		orders.run(Results.fromXML(new File(args[0])));
		System.out.println(orders);
	}

	private Map<Integer, Integer> nFoldsByOrder = new HashMap<Integer, Integer>();

	// just always call this "fold" in the code
	private Grouping normalizer = Grouping.fold();

	private double orderCutoff = 0.5;

	private double tmScoreCutoff = 0.4;

	public void run(Results census) {

		OrderDecider decider = new OrderDecider();
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);

		Map<String, Integer> nDomainsInFolds = new HashMap<String, Integer>();

		for (Result result : census.getData()) {
			try {
				ScopDomain domain = scop.getDomainByScopID(result.getScopId());
				if (domain == null) {
					logger.error(result.getScopId() + " is null");
				}
				result.setClassification(domain.getClassificationId());
				String fold = normalizer.group(result);
				StatUtils.plus(nDomainsInFolds, fold);
				decider.add(result);
			} catch (RuntimeException e) {
				logger.warn("Failed on " + result.getScopId(), e);
			}
		}

		for (String fold : nDomainsInFolds.keySet()) {
			int order = decider.getConsensusOrder(fold);
			StatUtils.plus(nFoldsByOrder, order);
		}

	}

	public Map<Integer, Integer> getnFoldsByOrder() {
		return nFoldsByOrder;
	}

	public void setNormalizer(Grouping normalizer) {
		this.normalizer = normalizer;
	}

	public void setOrderCutoff(double orderCutoff) {
		this.orderCutoff = orderCutoff;
	}

	public void setTmScoreCutoff(double tmScoreCutoff) {
		this.tmScoreCutoff = tmScoreCutoff;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("order\tN " + normalizer + StatUtils.NEWLINE);
		for (Map.Entry<Integer, Integer> entry : nFoldsByOrder.entrySet()) {
			sb.append(entry.getKey() + "\t" + entry.getValue() + StatUtils.NEWLINE);
		}
		return sb.toString();
	}

}
