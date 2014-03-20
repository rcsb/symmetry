package org.biojava3.structure.align.symm.census3.analysis.ligands;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.ElementType;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census3.stats.StructureClassificationGrouping;
import org.biojava3.structure.align.symm.census3.stats.CensusStatUtils;

/**
 * Prints statistics about {@link LigandList LigandLists}.
 * 
 * @author dmyersturnbull
 */
public class LigandStats {

	private static final LigandMatcher HIGH_OX_METALLIC_MATCHER = LigandMatcher.and(
			LigandMatcher.hasOxidationMagnitude(6), LigandMatcher.hasMetal());

	private static final LigandMatcher INTERESTING = LigandMatcher.and(LigandMatcher.hasPeriod(3), LigandMatcher
			.hasPeriod(8).not(), LigandMatcher.or(LigandMatcher.hasElementType(ElementType.TRANSITION_METAL),
					LigandMatcher.hasElementType(ElementType.POST_TRANSITION_METAL),
					LigandMatcher.hasElementType(ElementType.METALLOID)));

	private static final Logger logger = LogManager.getLogger(LigandStats.class.getName());

	private static final LigandMatcher METALLIC_NON_SALT = LigandMatcher.and(LigandMatcher.hasMetal(), LigandMatcher
			.isTrueSalt().not());

	private static String formatP(double x) {
		return formatD(x * 100) + "%";
	}

	private static String formatD(double x) {
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(0);
		return nf.format(x);
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + LigandFinder.class.getSimpleName() + " ligand-file.xml");
			return;
		}
		LigandList list = LigandList.fromXml(new File(args[0]));

		System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		System.out.println("+++++++++++++++++++++++NEAR CENTROID+++++++++++++++++++++++");
		System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		for (int i = 10; i > 0; i--) {
			LigandStats.printNormalizedStats(list, i, Double.POSITIVE_INFINITY);
		}
		System.out.println();

		System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		System.out.println("++++++++++++++++++++++ALONG INTERFACE++++++++++++++++++++++");
		System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		for (int i = 10; i > 0; i--) {
			LigandStats.printNormalizedStats(list, Double.POSITIVE_INFINITY, i);
		}

	}

	private static void printLine(String description, DescriptiveStatistics stats) {
		System.out.println(String.format("%-30s", description + ": ") + String.format("%8.1f", stats.getMean()) + "\t(stddev=" + String.format("%3.1f", stats.getStandardDeviation()) + ")");
	}
	
	public static void printNormalizedStats(LigandList list, double maxRadius, double maxAxisDistance) {
		LigandStats stats = new LigandStats();
		stats.setMaxRadius(maxRadius);
		stats.setMaxAxisDistance(maxAxisDistance);
		System.out.println("===============radius=" + maxRadius + "===============distance_to_axis=" + maxAxisDistance + "===============");
		System.out.println(String.format("%-30s", "!total superfamilies: ") + String.format("%8d", stats.getNSuperfamiliesTotal(list)));
		DescriptiveStatistics totalStats = stats.getNormalizedFractionMatching(list, LigandMatcher.everything());
		printLine("within constraints", totalStats);
		printLine("contains iron", stats.getNormalizedFractionMatching(list, LigandMatcher.hasElement(Element.Fe)));
		printLine("contains metal", stats.getNormalizedFractionMatching(list, LigandMatcher.hasMetal()));
		printLine("total |oxidation state| > 6",  stats.getNormalizedFractionMatching(list, LigandMatcher.hasTotalOxidationMagnitude(6)));
		printLine("any |oxidation state| > 6",  stats.getNormalizedFractionMatching(list, LigandMatcher.hasOxidationMagnitude(6)));
		printLine("contains non-salt metal",  stats.getNormalizedFractionMatching(list, METALLIC_NON_SALT));
		printLine("high-oxidation metallic",  stats.getNormalizedFractionMatching(list, HIGH_OX_METALLIC_MATCHER));
		printLine("is interesting",  stats.getNormalizedFractionMatching(list, INTERESTING));
	}

	/**
	 * @deprecated
	 */
	@Deprecated
	public static void printPerLigandStats(LigandList list, double maxRadius, double maxAxisDistance) {
		LigandStats stats = new LigandStats();
		stats.setMaxRadius(maxRadius);
		stats.setMaxAxisDistance(maxAxisDistance);
		System.out.println("===========================radius=" + maxRadius + "=================================");
		System.out.println("All: " + stats.getNLigandsMatching(list, LigandMatcher.everything()));
		System.out.println("Metallic: " + stats.getNLigandsMatching(list, LigandMatcher.hasMetal()));
		System.out.println("Metal: " + stats.getNLigandsMatching(list, LigandMatcher.isOnlyMetal()));
		System.out.println("Salt: " + stats.getNLigandsMatching(list, LigandMatcher.isTrueSalt()));
		System.out.println("Organic metal: " + stats.getNLigandsMatching(list, LigandMatcher.isOrganicMetal()));
		System.out.println("Has iron: " + stats.getNLigandsMatching(list, LigandMatcher.hasElement(Element.Fe)));
		System.out
		.println("High oxidation: " + stats.getNLigandsMatching(list, LigandMatcher.hasOxidationMagnitude(6)));
		System.out.println("High oxidation metallic: " + stats.getNLigandsMatching(list, HIGH_OX_METALLIC_MATCHER));
	}

	private double maxRadius = 5;

	private double maxAxisDistance = Double.POSITIVE_INFINITY;
	
	private int robustnessIterations = 200;

	private StructureClassificationGrouping normalizer = StructureClassificationGrouping.superfamily();

	/**
	 * @deprecated
	 */
	@Deprecated
	public int getNDomainsMatching(LigandList list, LigandMatcher matcher) {
		int n = 0;
		for (LigandsOfStructure inStruct : list.getData().values()) {
			for (CensusLigand ligand : inStruct.getLigands()) {
				if (ligand.getDistanceToCentroid() <= maxRadius) {
					if (matcher.matches(ligand)) {
						n++;
						break;
					}
				}
			}
		}
		return n;
	}

	/**
	 * @deprecated
	 */
	@Deprecated
	public int getNLigandsMatching(LigandList list, LigandMatcher matcher) {
		int n = 0;
		for (LigandsOfStructure inStruct : list.getData().values()) {
			for (CensusLigand ligand : inStruct.getLigands()) {
				if (ligand.getDistanceToCentroid() <= maxRadius) {
					if (matcher.matches(ligand)) {
						n++;
					}
				}
			}
		}
		return n;
	}

	public int getNSuperfamiliesTotal(LigandList list) {
		ScopDatabase scop = ScopFactory.getSCOP();
		Set<String> sfs = new HashSet<String>();
		for (Map.Entry<String, LigandsOfStructure> entry : list.getData().entrySet()) {
			ScopDomain domain = scop.getDomainByScopID(entry.getKey());
			if (domain == null) {
				logger.warn("Couldn't find domain " + domain);
				continue;
			}
			if (entry.getValue().isEmpty()) {
				continue;
			}
			String sf;
			try {
				sf = normalizer.group(domain);
			} catch (RuntimeException e) {
				throw new RuntimeException("Failed to get superfamily from domain " + domain, e);
			}
			if (sfs.contains(sf)) continue; // don't double-count SFs
			sfs.add(sf);
		}
		return sfs.size();
	}

	public DescriptiveStatistics getNormalizedFractionMatching(LigandList ligandList, LigandMatcher matcher) {
		
		ScopDatabase scop = ScopFactory.getSCOP();
		
		DescriptiveStatistics stats = new DescriptiveStatistics();
		
		// shuffle the order of domains
		List<String> shuffledKeyList = new ArrayList<String>(ligandList.size());
//		for (String s : ligandList.getData().keySet()) shuffledKeyList.add(s);
		for (Map.Entry<String, LigandsOfStructure> entry : ligandList.getData().entrySet()) {
			if (!entry.getValue().isEmpty()) {
				shuffledKeyList.add(entry.getKey());
			}
		}
		
		// do this repeatedly with different orders of domains each time
		for (int iter = 0; iter < robustnessIterations; iter++) {
			Collections.shuffle(shuffledKeyList);
			
			Set<String> sfs = new HashSet<String>();
			
			Set<String> containingSfs = new HashSet<String>();
			
			// make sure to loop over shuffledKeyList
			for (String domainKey : shuffledKeyList) {
				
				// get the domain and superfamily
				LigandsOfStructure inStruct = ligandList.get(domainKey);
				ScopDomain domain = scop.getDomainByScopID(domainKey);
				if (domain == null) {
					logger.warn("Couldn't find domain " + domain);
					continue;
				}
				String sf;
				try {
					sf = normalizer.group(domain);
				} catch (RuntimeException e) {
					throw new RuntimeException("Failed to get superfamily from domain " + domain, e);
				}
				
				// don't double-count SFs
				if (sfs.contains(sf)) continue;
				
				/*
				 * Now increment containingSfs UP TO MULTIPLICITY.
				 * Do this if even one Ligand in this structure matches.
				 */
				for (CensusLigand ligand : inStruct.getLigands()) {
					if (ligand.getDistanceToCentroid() <= maxRadius && ligand.getDistanceToAxis() <= maxAxisDistance) {
						sfs.add(sf);
						if (matcher.matches(ligand)) {
							containingSfs.add(sf);
							break;
						}
					}
				}
			}
			
			// record this into stats
			stats.addValue(containingSfs.size());
			
		}
//		logger.info("Standard deviation for " + matcher.toString() + " is " + stats.getStandardDeviation());
		return stats;
	}

	public int getNTotal(LigandList list) {
		return getNLigandsMatching(list, LigandMatcher.everything());
	}

	public void printStats(LigandList list) {
		int n = 0, nMetallic = 0;
		for (LigandsOfStructure inStruct : list.getData().values()) {
			boolean foundOne = false;
			boolean foundMetal = false;
			for (CensusLigand ligand : inStruct.getLigands()) {
				if (ligand.getDistanceToCentroid() <= maxRadius) {
					if (ligand.isMetallic()) foundMetal = true;
					foundOne = true;
				}
			}
			if (foundOne) n++;
			if (foundMetal) nMetallic++;
		}
		System.out.println(nMetallic + " / " + n + "\t" + CensusStatUtils.formatP((double) nMetallic / n));
	}

	public void setMaxRadius(double maxRadius) {
		this.maxRadius = maxRadius;
	}

	public void setMaxAxisDistance(double maxAxisDistance) {
		this.maxAxisDistance = maxAxisDistance;
	}

	public void setNormalizer(StructureClassificationGrouping normalizer) {
		this.normalizer = normalizer;
	}

	public void setRobustnessIterations(int robustnessIterations) {
		this.robustnessIterations = robustnessIterations;
	}

}
