package org.biojava.nbio.structure.align.symm.census3.stats.order;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDescription;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.census3.CensusResult;
import org.biojava.nbio.structure.align.symm.census3.CensusResultList;
import org.biojava.nbio.structure.align.symm.census3.CensusSignificance;
import org.biojava.nbio.structure.align.symm.census3.CensusSignificanceFactory;
import org.biojava.nbio.structure.align.symm.census3.stats.CensusStatUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Tabulate symmetry by order.
 * 
 * @author dmyersturnbull
 */
public class SymmetryOrder {
	private final static Logger logger = LoggerFactory.getLogger(SymmetryOrder.class);


	public static enum ExampleType {
		FOLD, SUPERFAMILY, FAMILY, DOMAIN;
	}

	public static class OrderInfo {
		private Map<String, Integer> nDomainsInFolds; // indexed by folds ONLY
		private Map<String, Set<String>> superfamilies; // indexed by classification (folds and superfamilies)
		private Map<String, Set<String>> families; // indexed by classification (folds and superfamilies)
		private Map<String, Set<String>> domains; // indexed by classification (folds, superfamilies, and families)
		private int order;

		public OrderInfo(int order) {
			this.order = order;
			nDomainsInFolds = new HashMap<String,Integer>();
			superfamilies = new HashMap<String, Set<String>>();
			families = new HashMap<String, Set<String>>();
			domains = new HashMap<String, Set<String>>();
		}

		public int getOrder() {
			return order;
		}

		public SortedSet<String> getDomainSet() {
			SortedSet<String> theSet = new TreeSet<String>(getComparator());
			for (Set<String> set : domains.values()) theSet.addAll(set);
			return theSet;
		}

		public SortedSet<String> getFamilySet() {
			SortedSet<String> theSet = new TreeSet<String>(getComparator());
			for (Set<String> set : families.values()) theSet.addAll(set);
			return theSet;
		}

		public SortedSet<String> getSuperfamilySet() {
			SortedSet<String> theSet = new TreeSet<String>(getComparator());
			for (Set<String> set : superfamilies.values()) {
				theSet.addAll(set);
			}
			return theSet;
		}

		private Comparator<String> getComparator() {
		return new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					if (!nDomainsInFolds.containsKey(o1) || !nDomainsInFolds.containsKey(o2))
						return o1.compareTo(o2);
					return nDomainsInFolds.get(o2).compareTo(nDomainsInFolds.get(o1));
				}
			};
		}
		
		public SortedSet<String> getFoldSet() {
			SortedSet<String> set = new TreeSet<String>(getComparator());
			for (String s : nDomainsInFolds.keySet()) set.add(s);
			return set;
		}

		public SortedMap<String, Integer> getnDomainsInFold() {
			Comparator<String> comp = new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					if (!nDomainsInFolds.containsKey(o1) || !nDomainsInFolds.containsKey(o2))
						return o1.compareTo(o2);
					return nDomainsInFolds.get(o2).compareTo(nDomainsInFolds.get(o1));
				}
			};
			SortedMap<String, Integer> counts = new TreeMap<String, Integer>(comp);
			counts.putAll(this.nDomainsInFolds);
			return counts;
		}

		public void plusFoldInstance(ScopDomain domain, ScopDatabase scop) {

			final String scopId = domain.getScopId();

			// determine fold
			int foldId = domain.getFoldId();
			ScopDescription foldDesc = scop.getScopDescriptionBySunid(foldId);
			String fold = foldDesc.getClassificationId();

			// determine superfamily
			int sfId = domain.getSuperfamilyId();
			ScopDescription sfDesc = scop.getScopDescriptionBySunid(sfId);
			String sf = sfDesc.getClassificationId();

			// determine family
			int familyId = domain.getFamilyId();
			ScopDescription famDesc = scop.getScopDescriptionBySunid(familyId);
			String family = famDesc.getClassificationId();

			CensusStatUtils.plus(nDomainsInFolds, fold);

			// plus superfamilies
			CensusStatUtils.plusSet(superfamilies, fold, sf);

			// plus families
			CensusStatUtils.plusSet(families, fold, family);
			CensusStatUtils.plusSet(families, sf, family);

			// plus domains
			CensusStatUtils.plusSet(domains, fold, scopId);
			CensusStatUtils.plusSet(domains, sf, scopId);
			CensusStatUtils.plusSet(domains, family, scopId);
		}

		@Override
		public String toString() {
			return toString(5);
		}

		/**
		 * Prints information about folds that have the symmetry order. Lists up to {@code limit} folds, ordered by number of domains with that symmetry.
		 * @param limit The maximum number of folds to list.
		 * @return
		 */
		public String toString(int limit) {
			SortedMap<String,Integer> counts = getnDomainsInFold();
			StringBuilder sb = new StringBuilder();
			sb.append(getDomainSet().size() + " domains, " + getFamilySet().size() + " families, " + getSuperfamilySet().size() + " superfamilies, " + getFoldSet().size() + " folds for order=" + order + ":" + CensusStatUtils.NEWLINE);
			int i = 0;
			sb.append("fold\tN domains\tN SFs" + CensusStatUtils.NEWLINE);
			for (Map.Entry<String, Integer> entry : counts.entrySet()) {
				if (i < limit) {
					sb.append(entry.getKey() + "\t" + entry.getValue() + "\t" + superfamilies.get(entry.getKey()).size()
							+ CensusStatUtils.NEWLINE);
				}
				i++;
			}
			return sb.toString();
		}
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
		SymmetryOrder orders = new SymmetryOrder(CensusResultList.fromXML(new File(args[0])));
		System.out.println(orders);
		System.out.println(orders.toTable(ExampleType.SUPERFAMILY, 16, "\t", "" + CensusStatUtils.NEWLINE, "\t"));
	}

	private Map<Integer, OrderInfo> orderInfos = new TreeMap<Integer, OrderInfo>();

	public SymmetryOrder(CensusResultList census) {
		ScopDatabase scop = ScopFactory.getSCOP();
		CensusSignificance sig = CensusSignificanceFactory.forCeSymmOrd();
		for (CensusResult result : census.getEntries()) {
			if (result.getAlignment() == null || result.getAxis() == null) {
				logger.warn("Skipping " + result.getId());
			}
			if (!sig.isSignificant(result))
				continue;
			int order;
			if (result.getOrder() > 1) {
				order = result.getOrder();
			} else {
				continue;
				//				order = result.getAxis().guessOrder();
			}
			ScopDomain domain = scop.getDomainByScopID(result.getId());
			if (domain == null) {
				logger.error(result.getId() + " is null");
			}
			if (!orderInfos.containsKey(order))
				orderInfos.put(order, new OrderInfo(order));
			orderInfos.get(order).plusFoldInstance(domain, scop);
		}
	}

	public String toTable(ExampleType exampleType, int exampleLimit, String tab, String newline, String comma) {
		StringBuilder sb = new StringBuilder();
		sb.append("order" + tab + "N folds" + tab + "N superfamilies" + tab + "N families" + tab + "N domains" + tab + "examples" + newline);
		for (OrderInfo info : orderInfos.values()) {
			sb.append(info.getOrder() + tab + info.getFoldSet().size() + tab + info.getSuperfamilySet().size() + tab + info.getFamilySet().size() + tab + info.getDomainSet().size() + tab);
			int i = 0;
			Set<String> examples = null;
			switch(exampleType) {
			case FOLD:
				examples = info.getFoldSet();
				break;
			case SUPERFAMILY:
				examples = info.getSuperfamilySet();
				break;
			case FAMILY:
				examples = info.getFamilySet();
				break;
			case DOMAIN:
				examples = info.getDomainSet();
				break;
			}
			for (String example : examples) {
				if (i >= exampleLimit) break;
				sb.append(example);
				i++;
				if (i < examples.size() && i < exampleLimit) sb.append(comma);
			}
			sb.append(newline);
		}
		return sb.toString();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (OrderInfo result : orderInfos.values()) {
			sb.append("====================== " + result.getOrder() + " ============================"
					+ CensusStatUtils.NEWLINE);
			sb.append(result.toString(10));
			sb.append("=====================================================" + CensusStatUtils.NEWLINE + CensusStatUtils.NEWLINE);
		}
		return sb.toString();
	}

}
