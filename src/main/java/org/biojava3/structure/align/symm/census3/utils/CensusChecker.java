package org.biojava3.structure.align.symm.census3.utils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusResultList;
import org.biojava3.structure.align.symm.census3.representatives.ScopSupport;

/**
 * Finds missing and extra SCOP Ids in a census file.
 * @author dmyersturnbull
 */
public class CensusChecker {

	private static final Logger logger = LogManager.getLogger(CensusChecker.class.getSimpleName());

	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + CensusChecker.class.getSimpleName() + " census-file.xml");
			return;
		}
		List<ScopDomain> list = new ArrayList<ScopDomain>();
		ScopSupport.getInstance().getAllDomainsUnder(ScopSupport.TRUE_SCOP_CLASSES, list, true);
		ScopSupport.nullifyInstance();
		Failures failures = findFailures(list, CensusResultList.fromXML(new File(args[0])));
		System.out.println(failures.toLineList());
	}
	
	public static class Failures {
		private List<String> missing = new ArrayList<String>();
		private List<String> extra = new ArrayList<String>();
		public List<String> getMissing() {
			return missing;
		}
		public List<String> getExtra() {
			return extra;
		}
		public void addMissing(String scopId) {
			missing.add(scopId);
		}
		public void addExtra(String scopId) {
			extra.add(scopId);
		}
		public String toLineList() {
			final String newline = System.getProperty("line.separator");
			StringBuilder sb = new StringBuilder();
			sb.append("---------------------------------------------------------------" + newline);
			sb.append("----------------------------MISSING----------------------------" + newline);
			sb.append("---------------------------------------------------------------" + newline);
			for (String scopId : getMissing()) sb.append(scopId + newline);
			sb.append("---------------------------------------------------------------" + newline);
			sb.append("-----------------------------EXTRA-----------------------------" + newline);
			sb.append("---------------------------------------------------------------" + newline);
			for (String scopId : getExtra()) sb.append(scopId + newline);
			return sb.toString();
		}
		@Override
		public String toString() {
			return missing.size() + " missing; " + extra.size() + " extra";
		}
	}
	
	public static Failures findFailures(List<ScopDomain> domains, CensusResultList census) {
		
		Failures failures = new Failures();
		
		// domains in list
		logger.info("Generating list of should-have domains...");
		List<String> scopIdsShouldHave = new ArrayList<String>(domains.size());
		for (ScopDomain domain : domains) scopIdsShouldHave.add(domain.getScopId());
		
		// domains in census
		logger.info("Generating list of domains in the census...");
		List<String> scopIdsHave = new ArrayList<String>(census.size());
		for (CensusResult result : census.getEntries()) scopIdsHave.add(result.getId());
		
		// count extra
		logger.info("Finding extra domains...");
		for (String scopId : scopIdsHave) {
			if (!scopIdsShouldHave.contains(scopId)) {
				failures.addExtra(scopId);
			}
		}
		
		// count missing
		logger.info("Finding missing domains...");
		for (String scopId : scopIdsShouldHave) {
			if (!scopIdsHave.contains(scopId)) {
				failures.addMissing(scopId);
			}
		}
		
		return failures;
		
	}

}
