package org.biojava3.structure.align.symm.census3.stats.order;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusResultList;
import org.biojava3.structure.align.symm.census3.stats.StructureClassificationGrouping;
import org.biojava3.structure.align.symm.census3.stats.CensusStatUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Tabulate symmetry by order, where symmetry order is determined by consensus over domains in a SCOP category.
 * 
 * @author dmyersturnbull
 */
public class FoldOrder {
	private final static Logger logger = LoggerFactory.getLogger(FoldOrder.class);

	public void setConsensusDecider(ConsensusDecider consensusDecider) {
		this.consensusDecider = consensusDecider;
	}

	private ConsensusDecider consensusDecider;
	
	public FoldOrder() {
		consensusDecider = ErrorMatrixDecider.fromMatrixFile();
	}

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
			orders.setNormalizer(StructureClassificationGrouping.byName(args[1]));
		}
		orders.run(CensusResultList.fromXML(new File(args[0])));
		System.out.println(orders);
	}

	private Map<Integer, Integer> nFoldsByOrder = new HashMap<Integer, Integer>();

	// just always call this "fold" in the code
	private StructureClassificationGrouping normalizer = StructureClassificationGrouping.fold();
	
	private StructureClassificationGrouping exampler = StructureClassificationGrouping.superfamily();

	public void setExampler(StructureClassificationGrouping exampler) {
		this.exampler = exampler;
	}

	private double tmScoreCutoff = 0.4;

	private Map<String, Map<String,Integer>> examplesByFold = new HashMap<String, Map<String,Integer>>();
	
	public void run(CensusResultList census) {

		OrderHelper decider = new OrderHelper(normalizer, tmScoreCutoff, consensusDecider);
		ScopDatabase scop = ScopFactory.getSCOP();

		Map<String, Integer> nDomainsInFolds = new HashMap<String, Integer>();

		for (CensusResult result : census.getEntries()) {
			try {
				ScopDomain domain = scop.getDomainByScopID(result.getId());
				if (domain == null) {
					logger.error(result.getId() + " is null");
				}
				String fold = normalizer.group(result);
				String superfamily = exampler.group(result);
				CensusStatUtils.plus(nDomainsInFolds, fold);
				if (!examplesByFold.containsKey(fold)) {
					examplesByFold.put(fold, new HashMap<String, Integer>());
				}
				CensusStatUtils.plus(examplesByFold.get(fold), superfamily);
				decider.add(result);
			} catch (RuntimeException e) {
				logger.warn("Failed on " + result.getId(), e);
			}
		}

		for (String fold : nDomainsInFolds.keySet()) {
			int order = decider.getConsensusOrder(fold);
			CensusStatUtils.plus(nFoldsByOrder, order);
		}

	}

	public Map<Integer, Integer> getnFoldsByOrder() {
		return nFoldsByOrder;
	}

	public void setNormalizer(StructureClassificationGrouping normalizer) {
		this.normalizer = normalizer;
	}

	public void setTmScoreCutoff(double tmScoreCutoff) {
		this.tmScoreCutoff = tmScoreCutoff;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("order\tN " + normalizer + CensusStatUtils.NEWLINE);
		for (Map.Entry<Integer, Integer> entry : nFoldsByOrder.entrySet()) {
			sb.append(entry.getKey() + "\t" + entry.getValue() + CensusStatUtils.NEWLINE);
		}
		return sb.toString();
	}

}
