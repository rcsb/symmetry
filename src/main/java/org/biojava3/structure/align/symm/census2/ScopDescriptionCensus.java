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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.representatives.ScopSupport;

/**
 * A census that takes a list of sun Ids.
 * @author dmyerstu
 */
public class ScopDescriptionCensus extends Census {

	private static final Logger logger = LogManager.getLogger(ScopDescriptionCensus.class.getPackage().getName());

	protected int[] sunIds;

	public static void buildDefault(File censusFile, int[] sunIds) {
		try {
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			ScopDescriptionCensus census = new ScopDescriptionCensus(maxThreads, sunIds);
			census.setOutputWriter(censusFile);
			census.setCache(new AtomCache());
			census.run();
			System.out.println(census);
		} catch (RuntimeException e) {
			logger.fatal(e.getMessage(), e);
		}
	}

	public static void main(String[] args) {
		if (args.length < 2) {
			System.err.println("Usage: " + ScopDescriptionCensus.class.getSimpleName() + " output-census-file sun-id-1 [sun-id-2 sun-id-3 ...]");
			return;
		}
		ScopFactory.setScopDatabase(ScopFactory.getSCOP(ScopFactory.VERSION_1_75A));
		final File censusFile = new File(args[0]);
		int[] sunIds = new int[args.length - 1];
		StringBuilder sb = new StringBuilder();
		for (int i = 1; i < args.length; i++) {
			Integer sunId = ScopSupport.getInstance().getSunId(args[i]);
			if (sunId == null) throw new IllegalArgumentException("Couldn't find " + args[i]);
			sb.append(sunId + ",");
			sunIds[i - 1] = sunId;
		}
		logger.info("Running on " + sb.toString().substring(0, sb.toString().length()-1));
		buildDefault(censusFile, sunIds);
	}

	public ScopDescriptionCensus(int maxThreads, int[] sunIds) {
		super(maxThreads);
		this.sunIds = sunIds;
	}

	@Override
	protected List<ScopDomain> getDomains() {
		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		ScopSupport.getInstance().getAllDomainsUnder(sunIds, domains);
		return domains;
	}

}
