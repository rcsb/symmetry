package org.biojava3.structure.align.symm.census2.analysis;

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
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

/**
 * Classify symmetry by Enzyme Commission number.
 * @author dmyerstu
 */
public class ECFinder {

	private static final Logger logger = LogManager.getLogger(ECFinder.class.getName());

	private Set<String> ecsByUnknownDomain = new HashSet<String>();
	private Map<String,String> ecsBySymmDomain = new HashMap<String,String>();
	private Map<String,String> ecsByAsymmDomain = new HashMap<String,String>();

	public static ECFinder fromTabbed(File file) throws IOException {
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(file));
			return fromTabbed(br);
		} finally {
			if (br != null) br.close();
		}
	}

	public static ECFinder fromTabbed(BufferedReader br) throws IOException {
		ECFinder corr = new ECFinder();
		String line = "";
		boolean onSymm = true;
		while ((line = br.readLine()) != null) {
			if (line.startsWith("<")) return null; // means XML
			if (line.startsWith("=")) {
				onSymm = false;
				continue;
			}
			String[] parts = line.split("\t");
			if (parts[1] != "-") {
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
		rebuild(Results.fromXML(file), output);
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
	
	public void rebuild(Results census, File output) {

		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
		Significance sig = SignificanceFactory.rotationallySymmetricSmart();

		int i = 0;
		for (Result result : census.getData()) {

			try {

				String scopId = result.getScopId();

				if (this.ecsByAsymmDomain.containsKey(scopId) || this.ecsBySymmDomain.containsKey(scopId) || ecsByUnknownDomain.contains(scopId)) {
					continue;
				}

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
						ecsBySymmDomain.put(scopId, ec);
					} else {
						ecsByAsymmDomain.put(scopId, ec);
					}

				} else if (ecs.size() > 1) {
					logger.info("Found different EC numbers for " + domain.getScopId()); // technically, this doesn't mean anything's wrong
				} else {
					logger.debug("Didn't find EC for " + scopId);
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
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String,String> entry : ecsBySymmDomain.entrySet()) {
			sb.append(entry.getKey() + "\t" + entry.getValue() + StatUtils.NEWLINE);
		}
		sb.append("=============================================" + StatUtils.NEWLINE);
		for (Map.Entry<String,String> entry : ecsByAsymmDomain.entrySet()) {
			sb.append(entry.getKey() + "\t" + entry.getValue() + StatUtils.NEWLINE);
		}
		for (String scopId : ecsByUnknownDomain) {
			sb.append(scopId + "\t-" + StatUtils.NEWLINE);
		}
		return sb.toString();
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: " + ECFinder.class.getSimpleName() + " census-file.xml ecs-file.tsv");
			return;
		}
		ECFinder ecs;
		if (new File(args[1]).exists()) {
			ecs = ECFinder.fromTabbed(new File(args[1]));
		} else {
			ecs = new ECFinder();
		}
		if (new File(args[0]).exists()) {
			ecs.rebuild(new File(args[0]), new File(args[2]));
		}
		System.out.println("============List of EC numbers of domains============");
		System.out.println("=====================================================" + StatUtils.NEWLINE);
		System.out.println("===================EC numbers level 0================");
//		ecs.printComparison(0, 10);
		System.out.println("=====================================================" + StatUtils.NEWLINE);
		System.out.println("===================EC numbers level 1================");
//		ecs.printComparison(1, 10);
		System.out.println("=====================================================" + StatUtils.NEWLINE);
	}

}
