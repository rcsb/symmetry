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
package org.biojava3.structure.align.symm.census3.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.CeSymm;

/**
 * A census that takes a file containing a line-by-line list of SCOP domains.
 * @author dmyersturnbull
 */
public class NamesCensus extends Census {

	private static final Logger logger = LogManager.getLogger(NamesCensus.class.getSimpleName());

	private List<ScopDomain> domains;

	public static void buildDefault(File censusFile, File lineByLine, final boolean doRefine) {
		// Alignment algorithm to actually run
		AlgorithmGiver algorithm = new AlgorithmGiver() {
			@Override
			public StructureAlignment getAlgorithm() {
				CeSymm ce = new CeSymm();
				ce.setRefineResult(doRefine);
				return ce;
			}
		};
		buildDefault(censusFile, lineByLine, algorithm);
	}

	public static void buildDefault(File censusFile, File lineByLine, AlgorithmGiver algorithm) {
		try {
			int maxThreads = 1;
			NamesCensus census = new NamesCensus(maxThreads);
			census.setOutputWriter(censusFile);
			census.domains = readNames(lineByLine);
			census.setPrintFrequency(10);
			census.setAlgorithm(algorithm);
			census.setRecordAlignmentMapping(true);
			AtomCache cache = new AtomCache();
			cache.setFetchFileEvenIfObsolete(true);
			census.setCache(cache);
			census.run();
			System.out.println(census);
		} catch (RuntimeException e) {
			logger.fatal(e.getMessage(), e);
		}
	}

	public static List<ScopDomain> readNames(File lineByLine) {
		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(lineByLine));
			String line = "";
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty()) continue;
				if(line.matches("^([;#]|//).*")) {
					continue; //comment
				}
				Scanner lineScanner = new Scanner(line);
				String token = lineScanner.next();
				lineScanner.close();
				ScopDomain domain = ScopFactory.getSCOP().getDomainByScopID(token);
				if (domain == null) {
					logger.error("No SCOP domain with id " + line + " was found");
				} else {
					domains.add(domain);
				}
			}
			br.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return domains;
	}

	public static void main(String[] args) {
		if (args.length != 2 && args.length != 3) {
			System.err.println("Usage: " + NamesCensus.class.getSimpleName() + " output-census-file line-by-line-input-names-file [do-refinement]");
			return;
		}
		final File censusFile = new File(args[0]);
		final File lineByLine = new File(args[1]);
		boolean doRefine = false;
		if (args.length > 2) {
			if (args[2].toLowerCase().equals("true") || args[2].toLowerCase().equals("refine")) {
				doRefine = true;
			}
		}
		ScopFactory.setScopDatabase(ScopFactory.getSCOP(ScopFactory.VERSION_1_75A,true));
		buildDefault(censusFile, lineByLine, doRefine);
	}

	public NamesCensus(int maxThreads) {
		super(maxThreads);
	}

	public NamesCensus(int maxThreads, List<ScopDomain> domains) {
		super(maxThreads);
		this.domains = domains;
	}

	@Override
	protected List<ScopDomain> getDomains() {
		return domains;
	}

}
