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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

/**
 * A census that takes a file containing a line-by-line list of SCOP domains.
 * @author dmyerstu
 */
public class NamesCensus extends Census {

	private List<ScopDomain> domains;

	public static void buildDefault(String pdbDir, File censusFile, File lineByLine) {
		try {
			ScopDatabase scop = Census.setBerkeleyScop(pdbDir);
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			NamesCensus census = new NamesCensus(maxThreads);
			census.setOutputWriter(censusFile);
			census.domains = new ArrayList<ScopDomain>();
			try {
				BufferedReader br = new BufferedReader(new FileReader(lineByLine));
				String line = "";
				while ((line = br.readLine()) != null) {
					if (line.trim().isEmpty()) continue;
					ScopDomain domain = scop.getDomainByScopID(line);
					if (domain == null) {
						logger.error("No SCOP domain with id " + line + " was found");
					} else {
						census.domains.add(domain);
					}
				}
				br.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			census.setCache(new AtomCache(pdbDir, false));
			census.run();
			System.out.println(census);
		} catch (RuntimeException e) {
			logger.warn(e.getMessage(), e);
		}
	}

	public static void main(String[] args) {
		final String pdbDir = args[0];
		final File censusFile = new File(args[1]);
		final File lineByLine = new File(args[2]);
		buildDefault(pdbDir, censusFile, lineByLine);
	}

	public NamesCensus(int maxThreads) {
		super(maxThreads);
	}

	public NamesCensus(int maxThreads, List<ScopDomain> domains) {
		super(maxThreads);
		this.domains = domains;
	}

	@Override
	protected Significance getSignificance() {
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return true;
			}
			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				return Census.getDefaultSignificance().isSignificant(protodomain, order, angle, afpChain);
			}
			@Override
			public boolean isSignificant(Result result) {
				return Census.getDefaultSignificance().isSignificant(result);
			}
		};
	}
	
	@Override
	protected List<ScopDomain> getDomains() {
		return domains;
	}

}
