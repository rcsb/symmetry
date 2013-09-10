package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.rcsb.RCSBDescription;
import org.biojava.bio.structure.rcsb.RCSBDescriptionFactory;
import org.biojava.bio.structure.rcsb.RCSBPolymer;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

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
					if (sig.isSignificant(result)) {
						symmEcs.put(scopId, ecs.first());
					} else {
						asymmEcs.put(scopId, ecs.first());
					}
				} else if (ecs.size() > 1) {
					logger.info("Found different EC numbers for " + domain.getScopId()); // technically, this doesn't mean anything's wrong
				}

				if (i > 0 && i % 10000 == 0) logger.debug("Working on #" + i);
				
			} catch (RuntimeException e) {
				logger.error(e);
			} finally {
				i++;
			}

		}
	}

	public void printComparison() {
		Map<Integer,Integer> symm = new HashMap<Integer,Integer>();
		for (String ec : symmEcs.values()) {
			int value = Integer.parseInt(ec.substring(0, 1));
			StatUtils.plus(symm, value);
		}
		Map<Integer,Integer> asymm = new HashMap<Integer,Integer>();
		for (String ec : asymmEcs.values()) {
			int value = Integer.parseInt(ec.substring(0, 1));
			StatUtils.plus(asymm, value);
		}
		for (int i = 1; i <= 6; i++) {
			double fSymm = 0, fAsymm = 0;
			if (symm.containsKey(i)) fSymm = (double) symm.get(new Integer(i)) / (double) symmEcs.size();
			if (asymm.containsKey(i)) fAsymm = (double) asymm.get(new Integer(i)) / (double) asymmEcs.size();
			System.out.println(i + "\t" + fSymm + "\t" + fAsymm);
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
		System.out.println(ecs);
		ecs.printComparison();
	}

}
