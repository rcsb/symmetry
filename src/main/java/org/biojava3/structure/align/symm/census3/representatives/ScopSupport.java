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
 * Created on 2013-06-23
 */
package org.biojava3.structure.align.symm.census3.representatives;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopNode;

/**
 * Utilities working with SCOP. A memory-intensive class that stores information from SCOP to save time. Uses a
 * singleton pattern, and has a {@link #nullifyInstance()} method to free memory.
 * 
 * @author dmyersturnbull
 */
public class ScopSupport {

	public static final int[] ALL_SCOP_CLASSES = new int[] { 46456, 48724, 51349, 53931, 56572, 56835, 56992, 57942, 58117, 58231, 58788 };
	public static final int[] TRUE_SCOP_CLASSES = new int[] { 46456, 48724, 51349, 53931, 56572, 56835 };

	public static void main(String[] args) throws IOException {
		ScopFactory.setScopDatabase(ScopFactory.LATEST_VERSION);
		for (int sunId : ALL_SCOP_CLASSES) {
			List<ScopDomain> domains = new ArrayList<ScopDomain>();
			ScopSupport.getInstance().getAllDomainsUnder(sunId, domains, true);
			String classId = ScopFactory.getSCOP().getScopDescriptionBySunid(sunId).getClassificationId();
			File file = new File("src/main/resources/scop_domains/" + classId + ".list");
			if (!file.exists()) {
				PrintWriter pw = null;
				try {
					pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
					for (ScopDomain domain : domains) {
						pw.println(domain.getScopId());
					}
				} finally {
					if (pw != null) pw.close();
				}
			}
		}
	}

	private static ScopSupport instance;
	private Map<ScopCategory, SoftReference<List<ScopDescription>>> categories = new HashMap<ScopCategory, SoftReference<List<ScopDescription>>>();

	private static final Logger logger = LogManager.getLogger(ScopSupport.class.getPackage().getName());

	private HashMap<String, Integer> descriptionsToSunIds = new HashMap<String, Integer>();

	public static ScopSupport getInstance() {
		if (instance == null) instance = new ScopSupport();
		return instance;
	}

	public ScopDescription getByIndex(ScopCategory category, int index) {
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> list = null;
		if (!categories.containsKey(category) || categories.get(category).get() == null) {
			list = scop.getByCategory(category);
			categories.put(category, new SoftReference<List<ScopDescription>>(list));
		} else {
			list = categories.get(categories.get(category)).get();
		}
		return list.get(index);
	}

	/**
	 * Set {@link #getInstance()} to null to free memory.
	 */
	public static void nullifyInstance() {
		instance = null;
		logger.info("Nullified ScopSupport");
	}

	private ScopSupport() {
		logger.info("Creating new ScopSupport...");
		ScopDatabase scop = ScopFactory.getSCOP();
		for (ScopCategory category : ScopCategory.values()) {
			if (category == ScopCategory.Species || category == ScopCategory.Px) continue;
			for (ScopDescription description : scop.getByCategory(category)) {
				descriptionsToSunIds.put(description.getClassificationId(), description.getSunID());
			}
		}
		logger.info("Finished creating ScopSupport");
	}


	/**
	 * General-use method that does not attempt to fetch a diverse set.
	 * @param sunId
	 * @param domains
	 */
	public void getAllDomainsUnder(int[] sunIds, List<ScopDomain> domains, boolean includeAllProteins) {
		for (int sunId : sunIds) {
			getAllDomainsUnder(sunId, domains, includeAllProteins);
		}
	}

	/**
	 * General-use method that does not attempt to fetch a diverse set.
	 * @param sunId
	 * @param domains
	 */
	public void getAllDomainsUnder(int sunId, List<ScopDomain> domains, boolean includeAllProteins) {

		// first just try to read from a file
		boolean success = readDomainsUnder(sunId, domains, includeAllProteins);
		if (success) return;

		final ScopDatabase scop = ScopFactory.getSCOP();
		final ScopDescription description = scop.getScopDescriptionBySunid(sunId);

		if (description == null) {
			throw new IllegalArgumentException("Couldn't find id " + sunId);
		}

		if (description.getCategory().equals(ScopCategory.Px)) { // base case
			if (includeAllProteins) {
				domains.addAll(scop.getScopDomainsBySunid(sunId));
			} else {
				domains.add(scop.getScopDomainsBySunid(sunId).get(0));
			}
		} else { // recurse
			final ScopNode node = scop.getScopNode(sunId);
			for (int s : node.getChildren()) {
				getAllDomainsUnder(s, domains, includeAllProteins);
			}
		}
	}

	/**
	 * Returns a list of all domains (in SCOP hierarchy, {@code dm}) that are descendants of the SCOP node corresponding
	 * to the Sun Id {@code sunId}. Unless {@code includeAllProteins} is set to true, only the <em>first</em> {@code Sp}
	 * and {@code Px} is chosen. When {@code repsPerSf} is less than the number of domains in a superfamily, this method
	 * tries to spread the chosen domains over every family in the superfamily. Specifically, the difference between the
	 * number of domains selected will not differ by more than 1 for any two families.
	 * <em>Warning: because this method looks over the entire SCOP subtree under the given node, it is very slow.</em>
	 * 
	 * @param sunId
	 *            The Sun Id of the ancestor node
	 * @param domains
	 *            A list of SCOP domains that will be appended to
	 * @param repsPerSf
	 *            The maximum number of domains to include per superfamily ({@code sf})
	 * @param includeAllProteins
	 *            Whether to include every {@code Sp} and {@code Px}; otherwise, choose only the first
	 */
	public void getDomainsUnder(int sunId, List<ScopDomain> domains, Integer repsPerSf, boolean includeAllProteins) {

		if (repsPerSf == null) repsPerSf = Integer.MAX_VALUE;

		final ScopDatabase scop = ScopFactory.getSCOP();

		final ScopDescription description = scop.getScopDescriptionBySunid(sunId);
		final ScopNode node = scop.getScopNode(sunId);
		final ScopCategory category = description.getCategory();

		if (category == ScopCategory.Class || category == ScopCategory.Fold) {

			// we're on class or fold, so recurse
			for (int s : node.getChildren()) {
				getDomainsUnder(s, domains, repsPerSf, includeAllProteins);
			}

		} else if (category == ScopCategory.Superfamily) {

			/*
			 * We would prefer having one domain from each family, as best as possible.
			 * To do this, make a list of queues of domains per family
			 * Then crawl through in the correct order
			 */

			int totalDomains = 0;
			// first, make the list
			List<Queue<ScopDomain>> domainsInFamilies = new ArrayList<Queue<ScopDomain>>();
			for (int familyId : node.getChildren()) {
				final ScopNode familyNode = scop.getScopNode(familyId);
				Queue<ScopDomain> inThisFamily = new LinkedList<ScopDomain>();
				for (int domainId : familyNode.getChildren()) {
					// just add the first species and Px available
					ScopDomain thisDomain = scop.getScopDomainsBySunid(domainId).get(0);
					inThisFamily.add(thisDomain);
					totalDomains++;
				}
				domainsInFamilies.add(inThisFamily);
			}

			logger.trace("Found " + totalDomains + " total domains from" + domainsInFamilies.size() + " families.");

			putDomainsFromFamilies(repsPerSf, totalDomains, domains, domainsInFamilies, description);

		} else if (category == ScopCategory.Family) {

			List<Queue<ScopDomain>> domainsInFamilies = new ArrayList<Queue<ScopDomain>>();
			Queue<ScopDomain> inThisFamily = new LinkedList<ScopDomain>();
			for (int domainId : node.getChildren()) {
				if (includeAllProteins) {
					inThisFamily.addAll(scop.getScopDomainsBySunid(domainId));
				} else {
					// just add the first species and Px available
					ScopDomain thisDomain = scop.getScopDomainsBySunid(domainId).get(0);
					inThisFamily.add(thisDomain);
				}
			}
			domainsInFamilies.add(inThisFamily);

			putDomainsFromFamilies(repsPerSf, node.getChildren().size(), domains, domainsInFamilies, description);

		} else {

			// just add the first Px if there are multiple
			// if we have a Px, this will always be the expected anyway
			ScopDomain domain = scop.getScopDomainsBySunid(sunId).get(0);
			domains.add(domain);

		}

	}

	/**
	 * Calls {@link #getDomainsUnder(int, List, Integer)} for each element in {@link sunIds}.
	 * 
	 * @see #getDomainsUnder(int, List, Integer)
	 */
	public void getDomainsUnder(int[] sunIds, List<ScopDomain> domains, Integer repsPerSf, boolean includeAllProteins) {
		for (int sunId : sunIds) {
			getDomainsUnder(sunId, domains, repsPerSf, includeAllProteins);
		}
	}

	/**
	 * Gets a Sun Id from a SCOP classification Id such as {@code k.38.1.1}, {@code d.40.1}, or {@code f}.
	 * 
	 * @param s
	 *            A classification Id; can also be a Sun Id in which case the value will be returned as long as it
	 *            exists
	 * @return The Sun Id, or null if it does not exist
	 */
	public Integer getSunId(String s) {
		ScopDatabase scop = ScopFactory.getSCOP();
		Integer sunId = null;
		try {
			sunId = Integer.parseInt(s);
			if (scop.getScopDescriptionBySunid(sunId) == null) sunId = null; // require that the Sun Id exists
		} catch (NumberFormatException e) {
			sunId = descriptionsToSunIds.get(s);
		}
		return sunId;
	}

	public Set<Integer> getSunIds(String[] strings) {
		Set<Integer> set = new HashSet<Integer>();
		for (String s : strings) {
			set.add(getSunId(s));
		}
		return set;
	}

	private boolean readDomainsUnder(int sunId, List<ScopDomain> domains, boolean includeAllProteins) {
		String classId = ScopFactory.getSCOP().getScopDescriptionBySunid(sunId).getClassificationId();
		File file = new File("src/main/resources/scop_domains/" + classId + ".list");
		if (!file.exists()) return false;
		logger.debug("Attempting to read " + file.getPath());
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(file));
			String line = "";
			while ((line = br.readLine()) != null) {
				ScopDomain domain = ScopFactory.getSCOP().getDomainByScopID(line.trim());
				if (domain == null) {
					logger.warn("Didn't find " + line.trim() + " in SCOP " + ScopFactory.getSCOP().getScopVersion() + "; perhaps the version is different?");
				} else {
					domains.add(domain);
				}
			}
		} catch (IOException e) {
			logger.warn("Couldn't read " + file.getPath() + "; re-constructing");
			return false;
		} finally {
			if (br != null) try {
				br.close();
			} catch (IOException e) {
				logger.warn("Couldn't close reader to " + file, e);
			}
		}
		return true;
	}

	private void putDomainsFromFamilies(int repsPerSf, int totalDomains, List<ScopDomain> domains,
			List<Queue<ScopDomain>> domainsInFamilies, ScopDescription description) {

		// okay, now add domains in a good order
		int preferredFamily = 0;
		for (int i = 0; i < repsPerSf; i++) {
			Queue<ScopDomain> queue;
			int nTried = 0;
			do {
				queue = domainsInFamilies.get(preferredFamily);
				if (preferredFamily < domainsInFamilies.size() - 1) {
					preferredFamily++;
				} else {
					preferredFamily = 0;
				}
				nTried++;
				if (nTried >= domainsInFamilies.size()) {
					break; // there are none left
				}
				if (queue.isEmpty()) continue;
			} while (queue.isEmpty());
			if (i >= totalDomains) {
				break;
			}
			if (!queue.isEmpty()) {
				final ScopDomain domain = queue.poll();
				domains.add(domain);
			} else {
				logger.warn("Missing 1 domain from a family in superfamily " + description.getClassificationId());
			}
		}

	}

}
