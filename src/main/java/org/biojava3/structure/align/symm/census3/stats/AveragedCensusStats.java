/*
 * BioJava development code
 * 
 * This code may be freely distributed and modified under the terms of the GNU Lesser General Public Licence. This
 * should be distributed with the code. If you do not have a copy, see:
 * 
 * http://www.gnu.org/copyleft/lesser.html
 * 
 * Copyright for this code is held jointly by the individual authors. These should be listed in @author doc comments.
 * 
 * For more information on the BioJava project and its aims, or to join the biojava-l mailing list, visit the home page
 * at:
 * 
 * http://www.biojava.org/
 * 
 * Created on 2013-03-10
 */
package org.biojava3.structure.align.symm.census3.stats;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusResultList;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Calculates statistics. To decide whether, say, a family is symmetric, it requires that the mean TM-score over each
 * domain in the family is at least the TM-score threshold, and that the mean order is at least the order threshold. To
 * decide whether a superfamily is symmetric, it requires that the mean TM-score over every <em>family</em> in that
 * superfamily is at least the TM-score threshold, and that the mean order is at least the order threshold. Similarly,
 * the symmetry of a fold is decided on the basis of mean values of its constituent superfamilies.
 * 
 * @author dmyersturnbull
 * @see {@link MajorityVotedCensusStats}, which uses majority voting instead
 */
public class AveragedCensusStats {

	private static double DEFAULT_ORDER_THRESHOLD = 0.5;

	private static double DEFAULT_TM_SCORE_THRESHOLD = 0.4;

	private final static Logger logger = LoggerFactory.getLogger(AveragedCensusStats.class);

	Map<String, Integer> nDomainsInFolds = new TreeMap<String, Integer>();

	Map<String, Integer> nInFoldsTotal = new TreeMap<String, Integer>();
	Map<String, Integer> nSymmInFolds = new TreeMap<String, Integer>();

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 3) {
			System.err.println("Usage: SmartStats census-file [tm-score-threshold] [order-threshold]");
			return;
		}
		File census = new File(args[0]);
		double tmScoreTreshold = DEFAULT_TM_SCORE_THRESHOLD;
		double orderTreshold = DEFAULT_ORDER_THRESHOLD;
		if (args.length > 1) tmScoreTreshold = Double.parseDouble(args[1]);
		if (args.length > 2) orderTreshold = Double.parseDouble(args[2]);
		AveragedCensusStats stats = new AveragedCensusStats(census, tmScoreTreshold, orderTreshold);
		System.out.println(stats);
	}

	public AveragedCensusStats(File census, double tmScoreTreshold, double orderTreshold) throws IOException {
		this(CensusResultList.fromXML(census), tmScoreTreshold, orderTreshold);
	}

	public AveragedCensusStats(CensusResultList census, double tmScoreTreshold, double orderTreshold) {

		Map<String, Double> tmInSfs = new TreeMap<String, Double>();
		Map<String, Integer> nInSfs = new TreeMap<String, Integer>();
		Map<String, Integer> nWithOrderInSfs = new TreeMap<String, Integer>();

		for (CensusResult result : census.getEntries()) {
			String classification = ScopFactory.getSCOP().getDomainByScopID(result.getId()).getClassificationId();
			double tmScore = 0;
			int order = 0;
			final String clas, fold, sf;
			try {
				String[] parts = classification.split("\\.");
				clas = parts[0];
				fold = parts[0] + "." + parts[1];
				sf = parts[0] + "." + parts[1] + "." + parts[2];
			} catch (ArrayIndexOutOfBoundsException e) {
				logger.error("Bad classification: " + classification, e);
				continue;
			}
			CensusStatUtils.plus(nInSfs, sf);
			try {
				tmScore = CensusResultPropertyFactory.tmScore().getProperty(result);
			} catch (PropertyUndefinedException e) {
				// okay, just leave as 0
			}
			try {
				order = (byte) CensusResultPropertyFactory.hasAnyOrder().getProperty(result);
			} catch (PropertyUndefinedException e) {
				// okay, just leave as 0
			}
			CensusStatUtils.plus(nDomainsInFolds, fold);
			CensusStatUtils.plus(nDomainsInFolds, clas);
			CensusStatUtils.plus(nDomainsInFolds, "overall");
			CensusStatUtils.plusD(tmInSfs, sf, tmScore);
			CensusStatUtils.plus(nWithOrderInSfs, sf, order);
		}

		for (Map.Entry<String, Integer> entry : nInSfs.entrySet()) {
			final String sf = entry.getKey();

			final String fold;
			final String clas;
			try {
				String[] parts = sf.split("\\.");
				fold = parts[0] + "." + parts[1];
				clas = parts[0];
			} catch (ArrayIndexOutOfBoundsException e) {
				logger.error("Bad classification: " + sf, e);
				continue;
			}
			final int nInSf = entry.getValue();
			final double fracWithOrder = nWithOrderInSfs.get(sf).doubleValue() / nInSf;
			final double meanTmScore = tmInSfs.get(sf) / nInSf;
			if (fracWithOrder >= orderTreshold && meanTmScore >= tmScoreTreshold) {
				CensusStatUtils.plus(nSymmInFolds, fold);
				CensusStatUtils.plus(nSymmInFolds, clas);
				CensusStatUtils.plus(nSymmInFolds, "overall");
			}
			CensusStatUtils.plus(nInFoldsTotal, fold);
			CensusStatUtils.plus(nInFoldsTotal, clas);
			CensusStatUtils.plus(nInFoldsTotal, "overall");
		}

	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String, Integer> entry : nInFoldsTotal.entrySet()) {
			final String fold = entry.getKey();
			final int total = entry.getValue();
			double fractionSymm;
			if (nSymmInFolds.get(fold) != null) {
				fractionSymm = (double) nSymmInFolds.get(fold) / (double) total;
			} else {
				fractionSymm = 0;
			}
			if (nDomainsInFolds.get(fold) >= 10 && fractionSymm >= 0.3 || fold.length() < 3 || fold.equals("overall")) {
				sb.append(fold + "\t" + nInFoldsTotal.get(fold) + "\t" + nDomainsInFolds.get(fold) + "\t"
						+ CensusStatUtils.formatP(fractionSymm) + CensusStatUtils.NEWLINE);
			}
		}
		return sb.toString();
	}

}
