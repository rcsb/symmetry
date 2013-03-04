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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
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

	static final Logger logger = Logger.getLogger(RandomDomains.class.getPackage().getName());

	static {
		BasicConfigurator.configure();
		logger.setLevel(Level.DEBUG);
	}
	
	private List<ScopDomain> domains;

	public List<ScopDomain> getDomains() {
		return domains;
	}

	public static void main(String[] args) {
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR, args[0]);
		ScopDatabase scop = ScopFactory.getSCOP();
		if (!scop.getClass().getName().equals(BerkeleyScopInstallation.class.getName())) { // for efficiency
			ScopFactory.setScopDatabase(new BerkeleyScopInstallation());
		}
		printRandomNames(Integer.parseInt(args[1]), System.out);
	}

	/**
	 * Prints a line-by-line list of SCOP ids randomly selected from the distribution of SCOP superfamilies. Only includes SCOP classes a, b, c, d, e, and f.
	 * @param number
	 * @param ps
	 */
	public static void printRandomNames(int number, PrintStream ps) {
		RandomDomains names = new RandomDomains(number);
		for (ScopDomain domain : names.domains) {
			ps.println(domain.getScopId());
		}
	}

	/**
	 * A recursive method that fills {@code descs} with a list of {@link ScopDescription ScopDescriptions} of the specified {@link ScopCategory} that are underneath {@code sunId} in the SCOP tree.
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

	public RandomDomains(int number) {
		this(number, ScopCategory.Superfamily, Arrays
				.asList(new Integer[] { 46456, 48724, 51349, 53931, 56572, 56835 }));
	}

	/**
	 * Creates a list of domains randomly selected according to the distribution of {@code category}, only including domains underneath an element in {@link classSunIds} in the SCOP tree.
	 * @param number The number of domains to generate
	 * @param category
	 * @param classSunIds A list of sun ids, required to be above {@link category} in the SCOP tree (need not actually be classes)
	 */
	public RandomDomains(int number, ScopCategory category, List<Integer> classSunIds) {
		Random random = new Random();
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
			logger.debug("Chose " + category.name() + " " + chosenCategory.getClassificationId() + " from " + allOfCategory.size() + " choices");
			List<ScopDescription> matchingDomainDescriptions = new ArrayList<ScopDescription>();
			getUnder(chosenCategory.getSunID(), ScopCategory.Domain, matchingDomainDescriptions);
			int domainDescriptionChoice = random.nextInt(matchingDomainDescriptions.size());
			ScopDescription chosenDescription = matchingDomainDescriptions.get(domainDescriptionChoice);
			logger.debug("Chose domain description " + chosenDescription.getClassificationId() + " from " + matchingDomainDescriptions.size() + " choices");
			List<ScopDomain> matchingDomains = scop.getScopDomainsBySunid(chosenDescription.getSunID());
			int domainChoice = random.nextInt(matchingDomains.size());
			ScopDomain chosenDomain = matchingDomains.get(domainChoice);
			logger.debug("Chose domain " + chosenDomain.getScopId() + " from " + matchingDomains.size() + " choices");
			domains.add(chosenDomain);
		}
	}

}
