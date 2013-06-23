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
 * Created on 2013-04-10
 */
package org.biojava3.structure.align.symm.census2.representatives;

import java.util.ArrayList;
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
 * 
 * @author dmyerstu
 */
public class SuperfamilyRepresentatives extends Representatives {

	private static final Logger logger = LogManager.getLogger(SuperfamilyRepresentatives.class.getPackage().getName());

	private List<ScopDomain> domains;

	private Integer repsPerSf;

	private int[] sunIds;

	public static void getDomainsUnder(int sunId, List<ScopDomain> domains, Integer repsPerSf) {

		if (repsPerSf == null) repsPerSf = Integer.MAX_VALUE;

		final ScopDatabase scop = ScopFactory.getSCOP();

		final ScopDescription description = scop.getScopDescriptionBySunid(sunId);
		final ScopNode node = scop.getScopNode(sunId);
		final ScopCategory category = description.getCategory();
		
		if (category == ScopCategory.Class || category == ScopCategory.Fold) {

			// we're on class or fold, so recurse
			for (int s : node.getChildren()) {
				getDomainsUnder(s, domains, repsPerSf);
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
				// just add the first species and Px available
				ScopDomain thisDomain = scop.getScopDomainsBySunid(domainId).get(0);
				inThisFamily.add(thisDomain);
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
	
	private static void putDomainsFromFamilies(int repsPerSf, int totalDomains, List<ScopDomain> domains, List<Queue<ScopDomain>> domainsInFamilies, ScopDescription description) {

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

	public SuperfamilyRepresentatives() {
		this(null);
	}

	public SuperfamilyRepresentatives(Integer numReps) {
		this(numReps, new int[] { 46456, 48724, 51349, 53931, 56572, 56835 });
	}

	/**
	 * Creates a new SuperfamilyRepresentives with {@code numReps} domains included for each Sun Id listed in {@link sunIds}.
	 * If {@code numReps} is less than the number of domains in a Sun Id, tries to spread the chosen domains over all families in that Sun Id.
	 * If {@code numReps == null}, then all domains from each Sun Id are always included.
	 * @param numReps
	 * @param sunIds
	 */
	public SuperfamilyRepresentatives(Integer numReps, int[] sunIds) {
		repsPerSf = numReps;
		this.sunIds = sunIds;
		domains = new ArrayList<ScopDomain>();
		for (int sunId : sunIds) {
			getDomainsUnder(sunId, domains, repsPerSf);
		}
	}

	@Override
	public List<ScopDomain> getDomains() {
		return domains;
	}

	public Integer getRepsPerSf() {
		return repsPerSf;
	}

	public int[] getSunIds() {
		return sunIds;
	}

}
