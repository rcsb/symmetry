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
package org.biojava3.structure.align.symm.census2.representatives;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

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

	private static ScopSupport instance;

	private static final Logger logger = LogManager.getLogger(ScopSupport.class.getPackage().getName());

	private HashMap<String, Integer> descriptionsToSunIds = new HashMap<String, Integer>();

	public static ScopSupport getInstance() {
		if (instance == null) instance = new ScopSupport();
		return instance;
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
	 * Returns a list of all domains (in SCOP hierarchy, {@code dm}) that are descendants of the SCOP node corresponding
	 * to the Sun Id {@code sunId}. Unless {@code includeAllProteins} is set to true, only the <em>first</em> {@code Sp}
	 * and {@code Px} is chosen. When {@code repsPerSf} is less than the number of domains in a superfamily, this method
	 * tries to spread the chosen domains over every family in the superfamily. Specifically, the difference between the
	 * number of domains selected will not differ by more than 1 for any two families.
	 * <em>Warning: because this method looks over the entire SCOP subtree under the given node, it is very slow.</em>
	 * 
	 * @param sunId
	 *            The Sun Id of the ancenstor node
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
