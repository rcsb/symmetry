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
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.scop.ScopDomain;

/**
 * A {@link Representatives} that contains domains under a list of SCOP Sun Ids, only including a specified number of
 * domains per SCOP superfamily.
 * 
 * @author dmyerstu
 */
public class SuperfamilyRepresentatives extends Representatives {

	private static final Logger logger = LogManager.getLogger(SuperfamilyRepresentatives.class.getPackage().getName());

	private List<ScopDomain> domains;

	private Integer repsPerSf;

	private int[] sunIds;

	public SuperfamilyRepresentatives() {
		this(null);
	}

	public SuperfamilyRepresentatives(Integer numReps) {
		this(numReps, new int[] { 46456, 48724, 51349, 53931, 56572, 56835 }, false);
	}

	/**
	 * Creates a new SuperfamilyRepresentives with {@code numReps} domains included for each Sun Id listed in
	 * {@link sunIds}. If {@code numReps} is less than the number of domains in a Sun Id, tries to spread the chosen
	 * domains over all families in that Sun Id. If {@code numReps == null}, then all domains from each Sun Id are
	 * always included.
	 * 
	 * @param numReps
	 * @param sunIds
	 * @param includeAllProteins
	 *            Whether to include every {@code Sp} and {@code Px}; otherwise, choose only the first
	 */
	public SuperfamilyRepresentatives(Integer numReps, int[] sunIds, boolean includeAllProteins) {
		repsPerSf = numReps;
		this.sunIds = sunIds;
		domains = new ArrayList<ScopDomain>();
		ScopSupport.getInstance().getDomainsUnder(sunIds, domains, repsPerSf, includeAllProteins);
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
