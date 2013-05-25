/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2013-03-03
 *
 */
package org.biojava3.structure.align.symm.census2.benchmark;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopNode;

/**
 * A set of random SCOP domains generated according to the distribution of some {@link ScopCategory}. The most basic
 * constructor generates the domains from the distribution of SCOP superfamilies on SCOP classes a, b, c, d, e, and f.
 * 
 * @author dmyerstu
 */
public class RandomDomains {

	private static final Logger logger = LogManager.getLogger(RandomDomains.class.getName());

	private List<ScopDomain> domains;

	/**
	 * Generates some random Berkeley SCOP domains from the distribution of superfamilies.
	 * @param args
	 *            <ol>
	 *            <li>PDB cache path</li>
	 *            <li>The number of domains to generate</li>
	 *            </ol>
	 */
	public static void main(String[] args) {
		if (args.length != 2) {
			System.err.println("Usage: RandomDomains pdb-dir number-of-domains");
			return;
		}
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR, args[0]);
		ScopDatabase scop = ScopFactory.getSCOP();
		if (!scop.getClass().getName().equals(BerkeleyScopInstallation.class.getName())) { // for efficiency
			ScopFactory.setScopDatabase(new BerkeleyScopInstallation());
		}
		System.out.println(printRandomNames(Integer.parseInt(args[1])));
	}

	/**
	 * Prints a line-by-line list of SCOP ids randomly selected from the distribution of SCOP superfamilies. Only
	 * includes SCOP classes a, b, c, d, e, and f.
	 * 
	 * @param number
	 * @param ps
	 */
	public static String printRandomNamesSimple(int number) {
		StringBuilder sb = new StringBuilder();
		RandomDomains names = new RandomDomains(number);
		for (ScopDomain domain : names.domains) {
			sb.append(domain.getScopId() + "\t" + linkifySymm(domain.getScopId()));
		}
		return sb.toString();
	}

	/**
	 * Prints a line-by-line list of SCOP ids randomly selected from the distribution of SCOP superfamilies. Only includes SCOP classes a, b, c, d, e, and f.
	 * @param number
	 * @param ps
	 */
	public static String printRandomNames(int number) {
		RandomDomains names = new RandomDomains(number);
		AtomCache cache = new AtomCache();
		StringBuilder sb = new StringBuilder();
		for (ScopDomain domain : names.domains) {
			sb.append(domain.getScopId() + "\n");
		}
		for (ScopDomain domain : names.domains) {
			String description = "";
			try {
				description = cache.getStructureForDomain(domain).getPDBHeader().getDescription();
				description = description.split("\\|")[1].substring(1);
			} catch (Exception e) {
				logger.error("Couldn't get PDB header or structure for " + domain.getScopId(), e);
			}
			sb.append("=hyperlink(\"" + linkify(domain.getScopId()) + "\";\"" + domain.getScopId() + "\")" + "\t" + description + "\n");
		}
		return (sb.toString());
	}

	private static String linkify(String scopId) {
		return "http://source.rcsb.org/jfatcatserver/showSymmetry.jsp?name1=" + scopId;
	}

	/**
	 * A recursive method that fills {@code descs} with a list of {@link ScopDescription ScopDescriptions} of the
	 * specified {@link ScopCategory} that are underneath {@code sunId} in the SCOP tree.
	 * 
	 * @param sunId
	 * @param category
	 * @param descs
	 */
	private static void getUnder(int sunId, ScopCategory category, List<ScopDescription> descs) {

		final ScopDatabase scop = ScopFactory.getSCOP();
		final ScopDescription description = scop.getScopDescriptionBySunid(sunId);

		if (description.getCategory().equals(category)) { // base case
			descs.add(scop.getScopDescriptionBySunid(sunId));
		} else { // recurse
			final ScopNode node = scop.getScopNode(sunId);
			for (int s : node.getChildren()) {
				getUnder(s, category, descs);
			}
		}
	}

	private static String linkifySymm(String scopId) {
		return "<a href=\"http://source.rcsb.org/jfatcatserver/showSymmetry.jsp?name1=" + scopId + "\">show</a>";
	}

	public RandomDomains(int number) {
		this(number, ScopCategory.Superfamily, Arrays
				.asList(new Integer[] { 46456, 48724, 51349, 53931, 56572, 56835 }));
	}

	/**
	 * See {@link #RandomDomains(int, ScopCategory, List, Random)}.
	 * 
	 * @param number
	 * @param category
	 * @param classSunIds
	 */
	public RandomDomains(int number, ScopCategory category, List<Integer> classSunIds) {
		this(number, category, classSunIds, new Random());
	}

	/**
	 * Creates a list of domains randomly selected according to the distribution of {@code category}, only including
	 * domains underneath an element in {@link classSunIds} in the SCOP tree.
	 * 
	 * @param number
	 *            The number of domains to generate
	 * @param category
	 * @param classSunIds
	 *            A list of sun ids, required to be above {@link category} in the SCOP tree (need not actually be
	 *            classes)
	 * @param random
	 *            A random-number generator
	 */
	public RandomDomains(int number, ScopCategory category, List<Integer> classSunIds, Random random) {
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> allOfCategory = new ArrayList<ScopDescription>();
		for (int classSunId : classSunIds) {
			getUnder(classSunId, category, allOfCategory);
		}
		logger.info("Found " + allOfCategory.size() + " categories in " + classSunIds.size() + " classes");
		domains = new ArrayList<ScopDomain>(number);
		for (int i = 0; i < number; i++) {
			logger.debug("Working on " + i);
			int categoryChoice = random.nextInt(allOfCategory.size());
			ScopDescription chosenCategory = allOfCategory.get(categoryChoice);
			logger.debug("Chose " + category.name() + " " + chosenCategory.getClassificationId() + " from "
					+ allOfCategory.size() + " choices");
			List<ScopDescription> matchingDomainDescriptions = new ArrayList<ScopDescription>();
			getUnder(chosenCategory.getSunID(), ScopCategory.Domain, matchingDomainDescriptions);
			int domainDescriptionChoice = random.nextInt(matchingDomainDescriptions.size());
			ScopDescription chosenDescription = matchingDomainDescriptions.get(domainDescriptionChoice);
			logger.debug("Chose domain description " + chosenDescription.getClassificationId() + " from "
					+ matchingDomainDescriptions.size() + " choices");
			List<ScopDomain> matchingDomains = scop.getScopDomainsBySunid(chosenDescription.getSunID());
			int domainChoice = random.nextInt(matchingDomains.size());
			ScopDomain chosenDomain = matchingDomains.get(domainChoice);
			logger.debug("Chose domain " + chosenDomain.getScopId() + " from " + matchingDomains.size() + " choices");
			domains.add(chosenDomain);
		}
	}

	public List<ScopDomain> getDomains() {
		return domains;
	}

}
