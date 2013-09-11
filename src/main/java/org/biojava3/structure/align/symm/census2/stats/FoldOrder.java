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
 * Tabulate symmetry by order, normalizing by a SCOP category.
 * 
 * @author dmyerstu
 */
public class FoldOrder {

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

			String fold = normalizer.group(result);

			StatUtils.plus(counts, fold);
			
			if (result.getAlignment() != null) {
				StatUtils.plusF(tmScoreMap, fold, result.getAlignment().getTmScore());
			}

			int order = 1;
			if (result.getOrder() != null) {
				order = result.getOrder();
				if (order < 0) order = 1; // it would be unfair to decrease the score
				StatUtils.plus(hasOrderMap, fold, order > 1 ? 1 : 0);
			}

			if (ordersMap.get(fold) == null) {
				ordersMap.put(fold, new HashMap<Integer,Integer>(7));
			}
			if (result.getAlignment() != null && result.getAlignment().getTmScore() >= tmScoreCutoff && order > 1) {
				StatUtils.plus(ordersMap.get(fold), order);
			} else {
//				StatUtils.plus(ordersMap.get(fold), 1);
			}

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

			// otherwise, return the mode order
			int maxCount = 0;
			int maximizingOrder = 0;
			for (int i = 2; i <= 8; i++) {
				Integer count = ordersMap.get(fold).get(i);
				if (count == null) count = 0;
				if (count > maxCount) {
					maxCount = count;
					maximizingOrder = i;
				}
			}
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

	private double orderCutoff = 0.00;

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
		sb.append("order\tN-folds" + StatUtils.NEWLINE);
		for (Map.Entry<Integer, Integer> entry : nFoldsByOrder.entrySet()) {
			sb.append(entry.getKey() + "\t" + entry.getValue() + StatUtils.NEWLINE);
		}
		return sb.toString();
	}

}
