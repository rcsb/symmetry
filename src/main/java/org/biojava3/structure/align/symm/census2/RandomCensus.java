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
import java.util.Collections;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * A census that runs on some specified number of domains in each superfamily, optionally randomized.
 * @author dmyersturnbull
 * @deprecated
 */
@Deprecated
public class RandomCensus extends Census {

	private static final Logger logger = LogManager.getLogger(RandomCensus.class.getName());

	private final int domainsPerSf;
	private final boolean shuffle;

	public static void buildDefault(File censusFile, int domainsPerSf, boolean shuffle) {
		try {
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			RandomCensus census = new RandomCensus(maxThreads, domainsPerSf, shuffle);
			census.setOutputWriter(censusFile);
			AtomCache cache = new AtomCache();
			cache.setFetchFileEvenIfObsolete(true);
			census.setCache(cache);
			census.run();
			System.out.println(census);
		} catch (RuntimeException e) {
			logger.warn(e.getMessage(), e);
		}
	}

	public static void main(String[] args) {
		if (args.length < 2) {
			System.err.println("Usage: " + RandomCensus.class.getSimpleName() + " output-census-file n-domains-per-superfamily [\"shuffle\"]");
			return;
		}
		ScopFactory.setScopDatabase(ScopFactory.getSCOP(ScopFactory.VERSION_1_75A));
		final File censusFile = new File(args[0]);
		final int domainsPerSf = Integer.parseInt(args[1]);
		boolean shuffle = false;
		if (args.length > 2) {
			if (args[2].toLowerCase().equals("shuffle") || args[2].toLowerCase().equals("true")) shuffle = true;
		}
		buildDefault(censusFile, domainsPerSf, shuffle);
	}

	public RandomCensus(int maxThreads, int domainsPerSf, boolean shuffle) {
		super(maxThreads);
		this.domainsPerSf = domainsPerSf;
		this.shuffle = shuffle;
	}

	@Override
	protected List<ScopDomain> getDomains() {
		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> superfamilies = scop.getByCategory(ScopCategory.Superfamily);
		for (ScopDescription superfamily : superfamilies) {
			List<ScopDomain> inSf = scop.getScopDomainsBySunid(superfamily.getSunID());
			if (shuffle) Collections.shuffle(inSf);
			int i = 0;
			for (ScopDomain domain : inSf) {
				if (i >= domainsPerSf) break;
				domains.add(domain);
				i++;
			}
		}
		return domains;
	}

}
