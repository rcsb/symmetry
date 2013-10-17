package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.ElementType;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.stats.Grouping;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

/**
 * Prints statistics about {@link LigandList LigandLists}.
 * 
 * @author dmyerstu
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
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(0);
		return nf.format(x * 100) + "%";
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

		for (int i = 10; i > 0; i--) {
			LigandStats.printNormalizedStats(list, i);
		}

	}

	public static void printNormalizedStats(LigandList list, int radius) {
		LigandStats stats = new LigandStats();
		stats.setMinRadius(radius);
		System.out.println("===========================radius=" + radius + "=================================");
		System.out.println("Has iron: "
				+ formatP(stats.getNormalizedFractionMatching(list, LigandMatcher.hasElement(Element.Fe))));
		System.out.println("Metallic: " + formatP(stats.getNormalizedFractionMatching(list, LigandMatcher.hasMetal())));
		System.out.println("High-oxidation: "
				+ formatP(stats.getNormalizedFractionMatching(list, LigandMatcher.hasOxidationMagnitude(6))));
		System.out.println("High-oxidation single-element: "
				+ formatP(stats.getNormalizedFractionMatching(list, LigandMatcher.and(LigandMatcher
						.hasOxidationMagnitude(6), LigandMatcher.nElements(1), LigandMatcher.hasElectronegativity(3.0f)
						.not()))));
		System.out.println("Non-salt metallic: "
				+ formatP(stats.getNormalizedFractionMatching(list, METALLIC_NON_SALT)));
		System.out.println("Interesting: " + formatP(stats.getNormalizedFractionMatching(list, INTERESTING)));
		System.out.println("High-oxidation metallic: "
				+ formatP(stats.getNormalizedFractionMatching(list, HIGH_OX_METALLIC_MATCHER)));
	}

	public static void printPerLigandStats(LigandList list, int radius) {
		LigandStats stats = new LigandStats();
		stats.setMinRadius(radius);
		System.out.println("===========================radius=" + radius + "=================================");
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

	private int minRadius = 5;

	private Grouping normalizer = Grouping.superfamily();

	public int getNDomainsMatching(LigandList list, LigandMatcher matcher) {
		int n = 0;
		for (StructureLigands inStruct : list.getData().values()) {
			for (Ligand ligand : inStruct.getLigands()) {
				if (ligand.getDistance() <= minRadius) {
					if (matcher.matches(ligand)) {
						n++;
						break;
					}
				}
			}
		}
		return n;
	}

	public int getNLigandsMatching(LigandList list, LigandMatcher matcher) {
		int n = 0;
		for (StructureLigands inStruct : list.getData().values()) {
			for (Ligand ligand : inStruct.getLigands()) {
				if (ligand.getDistance() <= minRadius) {
					if (matcher.matches(ligand)) {
						n++;
					}
				}
			}
		}
		return n;
	}

	public double getNormalizedFractionMatching(LigandList list, LigandMatcher matcher) {
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
		Set<String> sfs = new HashSet<String>();
		Set<String> containingSfs = new HashSet<String>();
		for (Map.Entry<String, StructureLigands> entry : list.getData().entrySet()) {
			StructureLigands inStruct = entry.getValue();
			ScopDomain domain = scop.getDomainByScopID(entry.getKey());
			String sf = normalizer.group(domain);
			if (sfs.contains(sf)) continue; // don't double-count SFs
			for (Ligand ligand : inStruct.getLigands()) {
				if (ligand.getDistance() <= minRadius) {
					sfs.add(sf);
					if (matcher.matches(ligand)) {
						containingSfs.add(sf);
						break;
					}
				}
			}
		}
		return (double) containingSfs.size() / sfs.size();
	}

	public int getNTotal(LigandList list) {
		return getNLigandsMatching(list, LigandMatcher.everything());
	}

	public void printStats(LigandList list) {
		int n = 0, nMetallic = 0;
		for (StructureLigands inStruct : list.getData().values()) {
			boolean foundOne = false;
			boolean foundMetal = false;
			for (Ligand ligand : inStruct.getLigands()) {
				if (ligand.getDistance() <= minRadius) {
					if (ligand.isMetallic()) foundMetal = true;
					foundOne = true;
				}
			}
			if (foundOne) n++;
			if (foundMetal) nMetallic++;
		}
		System.out.println(nMetallic + " / " + n + "\t" + StatUtils.formatP((double) nMetallic / n));
	}

	public void setMinRadius(int minRadius) {
		this.minRadius = minRadius;
	}

}
