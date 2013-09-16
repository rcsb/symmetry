package org.biojava3.structure.align.symm.census2.stats.order;

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
import org.biojava3.structure.align.symm.census2.stats.Grouping;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

/**
 * Tabulate symmetry by order, where symmetry order is determined by consensus over domains in a SCOP category.
 * 
 * @author dmyerstu
 */
public class FoldOrder {

	public void setConsensusDecider(ConsensusDecider consensusDecider) {
		this.consensusDecider = consensusDecider;
	}

	private ConsensusDecider consensusDecider;
	
	public FoldOrder() {
		consensusDecider = ErrorKernelDecider.fromMatrixFile();
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
	
	private Grouping exampler = Grouping.superfamily();

	public void setExampler(Grouping exampler) {
		this.exampler = exampler;
	}

	private double orderCutoff = 0.5;

	private double tmScoreCutoff = 0.4;

	private Map<String, Map<String,Integer>> examplesByFold = new HashMap<String, Map<String,Integer>>();
	
	public void run(Results census) {

		OrderHelper decider = new OrderHelper(normalizer, tmScoreCutoff, orderCutoff, consensusDecider);
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
				String superfamily = exampler.group(result);
				StatUtils.plus(nDomainsInFolds, fold);
				if (!examplesByFold.containsKey(fold)) {
					examplesByFold.put(fold, new HashMap<String, Integer>());
				}
				StatUtils.plus(examplesByFold.get(fold), superfamily);
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
