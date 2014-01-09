package org.biojava3.structure.align.symm.census3.analysis.ec;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
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
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusResultList;
import org.biojava3.structure.align.symm.census3.CensusSignificance;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;

/**
 * Classify symmetry by Enzyme Commission number.
 * Reads from and writes to tab-delimited files of the form:
 * <pre>
 * scopId	EC_number
 * scopId	-
 * </pre>
 * A hyphen-minus means the EC class could not be found.
 * @author dmyersturnbull
 */
public class CensusEnzymeFinder {

	private static final Logger logger = LogManager.getLogger(CensusEnzymeFinder.class.getName());

	private Set<String> ecsByUnknownDomain = new HashSet<String>();
	private Map<String,String> ecsBySymmDomain = new HashMap<String,String>();
	private Map<String,String> ecsByAsymmDomain = new HashMap<String,String>();

	/**
	 * Returns a set of SCOP Ids for which an EC number was not found.
	 */
	public Set<String> getEcsByUnknownDomain() {
		return ecsByUnknownDomain;
	}

	/**
	 * Returns a map of SCOP Ids and their corresponding EC numbers, restricted to the case of symmetric domains.
	 */
	public Map<String, String> getEcsByAsymmDomain() {
		return ecsByAsymmDomain;
	}

	/**
	 * Returns a map of SCOP Ids and their corresponding EC numbers, restricted to the case of asymmetric domains.
	 */
	public Map<String, String> getEcsBySymmDomain() {
		return ecsBySymmDomain;
	}

	/**
	 * Creates a new ECFinder from a tab-delimited file.
	 */
	public static CensusEnzymeFinder fromTabbed(File file) throws IOException {
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(file));
			return fromTabbed(br);
		} finally {
			if (br != null) br.close();
		}
	}

	public static CensusEnzymeFinder fromTabbed(BufferedReader br) throws IOException {
		CensusEnzymeFinder corr = new CensusEnzymeFinder();
		String line = "";
		boolean onSymm = true;
		while ((line = br.readLine()) != null) {
			if (line.startsWith("<")) return null; // means XML
			if (line.startsWith("=")) {
				onSymm = false;
				continue;
			}
			String[] parts = line.split("\t");
			if (!parts[1].equals("-")) {
				if (onSymm) {
					corr.ecsBySymmDomain.put(parts[0], parts[1]);
				} else {
					corr.ecsByAsymmDomain.put(parts[0], parts[1]);
				}
			} else {
				corr.ecsByUnknownDomain.add(parts[0]);
			}
		}
		return corr;
	}

	public void rebuild(File file, File output) throws IOException {
		buildFromCensus(CensusResultList.fromXML(file), output);
	}

	private void print(File output) {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(output);
			pw.print(toString());
			pw.flush();
		} catch (IOException e) {
			logger.error("Couldn't print to " + output);
		} finally {
			if (pw != null) pw.close();
		}
	}

	/**
	 * Attempts to find an EC number for every domain in {@code census}, skipping those which are already listed in this ECFinder.
	 * Prints to {@code output} periodically.
	 */
	public void buildFromCensus(CensusResultList census, File output) {

		ScopDatabase scop = ScopFactory.getSCOP();
		CensusSignificance sig = CensusSignificanceFactory.forCeSymmOrd();

		int i = 0;
		for (CensusResult result : census.getEntries()) {

			try {

				String scopId = result.getId();

				if (this.ecsByAsymmDomain.containsKey(scopId) || this.ecsBySymmDomain.containsKey(scopId) || ecsByUnknownDomain.contains(scopId)) {
					continue;
				}

				ScopDomain domain = scop.getDomainByScopID(scopId);
				if (domain == null) {
					logger.error(result.getId() + " is null");
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
						ecsBySymmDomain.put(scopId, ec);
					} else {
						ecsByAsymmDomain.put(scopId, ec);
					}

				} else if (ecs.size() > 1) {
					logger.info("Found different EC numbers for " + domain.getScopId()); // technically, this doesn't mean anything's wrong
				} else {
//					logger.debug("Didn't find EC for " + scopId);
					ecsByUnknownDomain.add(scopId);
				}

				if (i > 0 && i % 100 == 0) {
					print(output);
					logger.debug("Working on #" + i);
				}

			} catch (RuntimeException e) {
				e.printStackTrace();
				logger.error(e);
			} finally {
				i++;
			}

		}

	}

	@Override
	public String toString() {
		String newline = System.getProperty("line.separator");
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String,String> entry : ecsBySymmDomain.entrySet()) {
			sb.append(entry.getKey() + "\t" + entry.getValue() + newline);
		}
		sb.append("=============================================" + newline);
		for (Map.Entry<String,String> entry : ecsByAsymmDomain.entrySet()) {
			sb.append(entry.getKey() + "\t" + entry.getValue() + newline);
		}
		for (String scopId : ecsByUnknownDomain) {
			sb.append(scopId + "\t-" + newline);
		}
		return sb.toString();
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: " + CensusEnzymeFinder.class.getSimpleName() + " census-file.xml ecs-file.tsv");
			return;
		}
		CensusEnzymeFinder ecs;
		if (new File(args[1]).exists()) {
			ecs = CensusEnzymeFinder.fromTabbed(new File(args[1]));
		} else {
			ecs = new CensusEnzymeFinder();
		}
		ecs.rebuild(new File(args[0]), new File(args[1]));
		CensusEnzymeStats.printComparisonUnnormalized(1, 10, ecs.ecsBySymmDomain, ecs.ecsBySymmDomain);
		CensusEnzymeStats.printComparisonUnnormalized(2, 10, ecs.ecsBySymmDomain, ecs.ecsBySymmDomain);
	}

}
