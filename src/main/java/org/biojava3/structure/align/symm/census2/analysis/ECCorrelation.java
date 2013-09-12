package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.rcsb.RCSBDescription;
import org.biojava.bio.structure.rcsb.RCSBDescriptionFactory;
import org.biojava.bio.structure.rcsb.RCSBPolymer;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

/**
 * Classify symmetry by Enzyme Classification.
 * @author dmyerstu
 */
public class ECCorrelation {

	private static final Logger logger = LogManager.getLogger(ECCorrelation.class.getName());

	private Map<String,String> symmEcs = new HashMap<String,String>();
	private Map<String,String> asymmEcs = new HashMap<String,String>();

	public ECCorrelation(Results census) {

		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
		Significance sig = SignificanceFactory.rotationallySymmetricSmart();

		int i = 0;
		for (Result result : census.getData()) {

			try {

				String scopId = result.getScopId();
				ScopDomain domain = scop.getDomainByScopID(scopId);
				if (domain == null) {
					logger.error(result.getScopId() + " is null");
					continue;
				}

				// got a result; what's its EC?

				// we need to find the correct polymers corresponding to the domain
				// note that this still isn't perfect, since we don't know what part of the polymer actually does the function
				List<RCSBPolymer> polymers = new ArrayList<RCSBPolymer>();
				Set<String> chains = domain.getChains();
				RCSBDescription desc = RCSBDescriptionFactory.get(domain.getPdbId());
				for (RCSBPolymer polymer : desc.getPolymers()) {
					for (Character chain : polymer.getChains()) {
						if (chains.contains(String.valueOf(chain))) {
							polymers.add(polymer);
							break;
						}
					}
				}

				// get the EC numbers
				// use a set because we don't want > 1 just because we have duplicates
				NavigableSet<String> ecs = new TreeSet<String>();
				for (RCSBPolymer polymer : polymers) {
					String ec = polymer.getEnzClass();
					if (ec != null) ecs.add(ec);
				}

				if (ecs.size() == 1) {

					String ec = ecs.first();

					if (sig.isSignificant(result)) {
						symmEcs.put(scopId, ec);
					} else {
						asymmEcs.put(scopId, ec);
					}

				} else if (ecs.size() > 1) {
					logger.info("Found different EC numbers for " + domain.getScopId()); // technically, this doesn't mean anything's wrong
				}

				if (i > 0 && i % 10000 == 0) logger.debug("Working on #" + i);

			} catch (RuntimeException e) {
				e.printStackTrace();
				logger.error(e);
			} finally {
				i++;
			}

		}
	}

	private String getLabel(String ec, int level) {
		String[] parts = ec.split("\\.");
		String label = parts[0];
		for (int i = 1; i <= level; i++) {
			if (i >= parts.length) return null; // can happen if the EC number isn't fully specified (in fact, this is common)
			label += "." + parts[i];
		}
		return label;
	}

	public void printComparison(int level, int nExamples) {

		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);

		if (level < 0 || level > 3) throw new IllegalArgumentException("Level must be between 0 and 3, inclusive");
		Set<String> labels = new LinkedHashSet<String>();

		final Map<String,Map<String,Integer>> foldsInEcs = new HashMap<String,Map<String,Integer>>();

		Map<String,Integer> symm = new HashMap<String,Integer>();
		for (Map.Entry<String,String> entry : symmEcs.entrySet()) {

			final String scopId = entry.getKey();
			final String ec = entry.getValue();

			String label = getLabel(ec, level);
			if (label == null) continue;

			// record the fold
			if (!foldsInEcs.containsKey(label)) foldsInEcs.put(label, new HashMap<String,Integer>());
			ScopDomain domain = scop.getDomainByScopID(scopId);
			ScopDescription desc = scop.getScopDescriptionBySunid(domain.getFoldId());
			String fold = desc.getName();
			StatUtils.plus(foldsInEcs.get(label), fold);

			StatUtils.plus(symm, label);
			labels.add(label);
		}

		Map<String,Integer> asymm = new HashMap<String,Integer>();
		for (String ec : asymmEcs.values()) {
			String label = getLabel(ec, level);
			if (label == null) continue;
			StatUtils.plus(asymm, label);
			labels.add(label);
		}

		for (String label : labels) {
			double fSymm = 0, fAsymm = 0;
			if (symm.containsKey(label)) fSymm = (double) symm.get(label);
			if (asymm.containsKey(label)) fAsymm = (double) asymm.get(label);
			System.out.print(label + "\t" + StatUtils.formatP(fSymm / (fSymm+fAsymm)) + "\t" + (fSymm+fAsymm));
			final Map<String,Integer> map = foldsInEcs.get(label);
			if (map != null) {
				Comparator<String> comp = new Comparator<String>() {
					@Override
					public int compare(String o1, String o2) {
						if (!map.containsKey(o1) || !map.containsKey(o2)) return 0;
						return map.get(o2).compareTo(map.get(o1));
					}
				};
				SortedMap<String,Integer> examples = new TreeMap<String,Integer>(comp);
				examples.putAll(map);
				int i = 0;
				for (Map.Entry<String,Integer> entry : examples.entrySet()) {
					System.out.print("\t" + entry.getKey());
					i++;
					if (i > nExamples) break;
				}
			}
		}
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String,String> entry : symmEcs.entrySet()) {
			System.out.println(entry.getKey() + "\t" + entry.getValue());
		}
		return sb.toString();
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + ECCorrelation.class.getSimpleName() + " census-file.xml");
			return;
		}
		ECCorrelation ecs = new ECCorrelation(Results.fromXML(new File(args[0])));
		System.out.println("============List of EC numbers of domains============");
		System.out.println(ecs);
		System.out.println("=====================================================" + StatUtils.NEWLINE);
		System.out.println("===================EC numbers level 0================");
		ecs.printComparison(0, 10);
		System.out.println("=====================================================" + StatUtils.NEWLINE);
		System.out.println("===================EC numbers level 1================");
		ecs.printComparison(1, 10);
		System.out.println("=====================================================" + StatUtils.NEWLINE);
	}

}
