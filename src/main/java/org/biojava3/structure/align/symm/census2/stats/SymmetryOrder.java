package org.biojava3.structure.align.symm.census2.stats;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

public class SymmetryOrder {

	private static final Logger logger = LogManager.getLogger(SymmetryOrder.class.getName());

	public static class ClassificationOrder {
		private Map<String,Integer> counts;
		private Map<String,Set<String>> sfs;
		private int order;
		private int totalCount;
		public ClassificationOrder(int order) {
			this.order = order;
			counts = new TreeMap<String,Integer>();
			sfs = new HashMap<String,Set<String>>();
		}
		public Map<String, Integer> getCounts() {
			return counts;
		}
		public int getOrder() {
			return order;
		}
		public int getTotalCount() {
			return totalCount;
		}
		public void plusFoldInstance(ScopDomain domain, ScopDatabase scop) {
			totalCount++;
			int foldId = domain.getFoldId();
			ScopDescription foldDesc = scop.getScopDescriptionBySunid(foldId);
			String fold = foldDesc.getClassificationId();
			int sfId = domain.getSuperfamilyId();
			ScopDescription sfDesc = scop.getScopDescriptionBySunid(sfId);
			String sf = sfDesc.getClassificationId();
			StatUtils.plus(counts, fold);
			if (!sfs.containsKey(fold)) sfs.put(fold, new HashSet<String>());
			sfs.get(fold).add(sf);
		}
		@Override
		public String toString() {
			return toString(5);
		}
		public String toString(int limit) {
			Comparator<String> comp = new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					if (!counts.containsKey(o1) || !counts.containsKey(o2)) return o1.compareTo(o2);
					return counts.get(o2).compareTo(counts.get(o1));
				}
			};
			SortedMap<String,Integer> counts = new TreeMap<String,Integer>(comp);
			counts.putAll(this.counts);
			StringBuilder sb = new StringBuilder();
			sb.append(totalCount + " results for order=" + order + ":" + StatUtils.NEWLINE);
			int i = 0;
			for (Map.Entry<String,Integer> entry : counts.entrySet()) {
				if (i < limit) {
					sb.append(entry.getKey() + "\t" + entry.getValue() + "\t" + sfs.get(entry.getKey()).size() + StatUtils.NEWLINE);
				}
				i++;
			}
			return sb.toString();
		}
	}

	private Map<Integer,ClassificationOrder> results = new TreeMap<Integer,ClassificationOrder>();

	public SymmetryOrder(Results census) {
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
		Significance sig = SignificanceFactory.rotationallySymmetricSmart();
		for (Result result : census.getData()) {
			if (result.getAlignment() == null || result.getAxis() == null) {
				logger.warn("Skipping " + result.getScopId());
			}
			if (!sig.isSignificant(result)) continue;
			int order;
			if (result.getOrder() != null && result.getOrder() > 1) {
				order = result.getOrder();
			} else {
				order = result.getAxis().guessOrder();
			}
			if (order < 2) {
				System.err.println(result.toString());
				System.out.println(sig.isSignificant(result));
				throw new RuntimeException("AAAHAHAHHHHHHH");
			}
			ScopDomain domain = scop.getDomainByScopID(result.getScopId());
			if (domain == null) {
				logger.error(result.getScopId() + " is null");
			}
			if (!results.containsKey(order)) results.put(order, new ClassificationOrder(order));
			results.get(order).plusFoldInstance(domain, scop);
		}
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (ClassificationOrder result : results.values()) {
			sb.append("====================== " + result.getOrder() + " ============================" + StatUtils.NEWLINE);
			sb.append(result.toString(10));
			sb.append("=====================================================" + StatUtils.NEWLINE + StatUtils.NEWLINE);
		}
		return sb.toString();
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + SymmetryOrder.class.getSimpleName() + " census-file.xml");
			return;
		}
		SymmetryOrder orders = new SymmetryOrder(Results.fromXML(new File(args[0])));
		System.out.println(orders);
	}

}
