package org.biojava3.structure.align.symm.census2.stats;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.jama.Matrix;
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

	public void setConsensusDecider(ConsensusDecider consensusDecider) {
		this.consensusDecider = consensusDecider;
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
	 */
	public static class ErrorKernelDecider implements ConsensusDecider {

		private Matrix kernel;

		public ErrorKernelDecider(Matrix kernel) {
			super();
			this.kernel = kernel;
		}

		@Override
		public int decide(Map<Integer,Integer> countsByOrders) {
			
			double[] flows = new double[] {0, 0, 0, 0, 0, 0, 0};
			
			for (int i = 2; i <= 8; i++) { // correct
				for (int j = 2; j <= 8; j++) { // putative
					
					// the number of states we're allowed to use
					int m = i / j;
					
					/*
					 *  use Chapman–Kolmogorov equation to find m-state transition kernel
					 *  We want the m-state transition because we're only concerned with flow from i to j,
					 *  which we know requires EXACTLY m steps
					 */
					Matrix mStateTransition = kernel;
					for (int k = 1; k < m; k++) mStateTransition = mStateTransition.times(kernel);
					
					// get the number of "j"s, our putative actual "i"s
					Integer count = countsByOrders.get(j);
					if (count == null) count = 0;
					
					/*
					 * We want to make "putative" flow into "correct"
					 * Ex: Matrix[4,2] = 0.1, the probability we got 2 instead of 4
					 * Then we're using mStateTransition[4, 2] = 0.1
					 * And thus we allow 0.1 * count to flow from 2 into 4
					 * So this is the correct indexing
					 */
					double flow = mStateTransition.get(i - 2, j - 2);
					flows[i-2] += count * flow;
					
				}
			}
			
			double maxFlow = 0;
			int maximizingOrder = 0;
			for (int i = 2; i <= 8; i++) {
				double flow = flows[i-2];
				if (flow > maxFlow) {
					maxFlow = flow;
					maximizingOrder = i;
				}
			}
			return maximizingOrder;
		}
	}

	private static Matrix KERNEL;
	static {
		BufferedReader br = null;
		try {
			try {
				br = new BufferedReader(new FileReader("src/main/resources/error_kernel.matrix"));
				KERNEL = Matrix.read(br);
			} finally {
				if (br != null) br.close();
			}
		} catch (IOException e) {
			e.printStackTrace(); // might be ok
		}
	}

	private ConsensusDecider consensusDecider = new ErrorKernelDecider(KERNEL);

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
