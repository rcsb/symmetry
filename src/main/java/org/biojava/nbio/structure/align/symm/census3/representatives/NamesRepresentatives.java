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
 * Created on 2013-04-10
 *
 */
package org.biojava.nbio.structure.align.symm.census3.representatives;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * Representatives from a simple list of SCOP Ids.
 * @author dmyersturnbull
 */
public class NamesRepresentatives extends Representatives {

	@Override
	public List<ScopDomain> getDomains() {
		return domains;
	}

	public NamesRepresentatives(String... names) {
		this(Arrays.asList(names));
	}
	
	public NamesRepresentatives(List<String> names) {
		domains = new ArrayList<ScopDomain>();
		for (String name : names) {
			domains.add(ScopFactory.getSCOP().getDomainByScopID(name));
		}
	}
	
	private List<ScopDomain> domains;
}
