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
 * Created on 2013-02-18
 *
 */
package org.biojava3.structure.align.symm.census2;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopNode;

/**
 * A census that takes a list of sun Ids
 * @author dmyerstu
 */
public class ScopDescriptionCensus extends Census {

	static final Logger logger = Logger.getLogger(CensusJob.class.getPackage().getName());
	protected int[] sunIds;

	static {
		BasicConfigurator.configure();
	}

	public static void buildDefault(String pdbDir, File censusFile, int[] sunIds) {
		try {
			Census.setBerkeleyScop(pdbDir);
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			ScopDescriptionCensus census = new ScopDescriptionCensus(maxThreads, sunIds);
			census.setOutputWriter(censusFile);
			census.setCache(new AtomCache(pdbDir, false));
			census.run();
			System.out.println(census);
		} catch (RuntimeException e) {
			logger.fatal(e.getMessage(), e);
		}
	}

	public static void main(String[] args) {
		final String pdbDir = args[0];
		final File censusFile = new File(args[1]);
		int[] sunIds = new int[args.length - 2];
		for (int i = 2; i < args.length; i++) {
			sunIds[i - 2] = Integer.parseInt(args[i]);
		}
		buildDefault(pdbDir, censusFile, sunIds);
	}

	public ScopDescriptionCensus(int maxThreads, int[] sunIds) {
		super(maxThreads);
		this.sunIds = sunIds;
	}

	@Override
	protected List<ScopDomain> getDomains() {
		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		ScopFactory.setScopDatabase(new BerkeleyScopInstallation()); // TODO Why the hell do I have to do this again here?
		ScopDatabase scop = ScopFactory.getSCOP();
		for (int sunId : sunIds) {
			domains.addAll(scop.getScopDomainsBySunid(sunId));
			getDomainsUnder(sunId, domains);
		}
		return domains;
	}

	public static void getDomainsUnder(int sunId, List<ScopDomain> domains) {
		
		final ScopDatabase scop = ScopFactory.getSCOP();
		final ScopDescription description = scop.getScopDescriptionBySunid(sunId);
		
		if (description.getCategory().equals(ScopCategory.Domain)) { // base case
			domains.addAll(scop.getScopDomainsBySunid(sunId));
		} else { // recurse
			final ScopNode node = scop.getScopNode(sunId);
			for (int s : node.getChildren()) {
				getDomainsUnder(s, domains);
			}
		}
	}
	
}
