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
package org.biojava3.structure.align.symm.census3.representatives;

import java.util.HashMap;
import java.util.List;

import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * A list of {@link ScopDomain ScopDomains} to be used as representatives.
 * This class also implements a singleton pattern, although Representative objects aren't required to be used this way.
 * @author dmyersturnbull
 */
public abstract class Representatives {

	private static volatile Representatives defaultNames;
	private static Object defaultNamesLock = new Object();

	private HashMap<String, ScopDescription> superfamilies;

	public static Representatives get() {
		if (defaultNames == null) setDefault();
		return defaultNames;
	}

	public static void set(Representatives representatives) {
		defaultNames = representatives;
	}

	public static void setDefault() {
		synchronized (defaultNamesLock) {
			defaultNames = new SuperfamilyRepresentatives();
		}
	}

	public abstract List<ScopDomain> getDomains();

	public HashMap<String, ScopDescription> getSuperfamilies() {
		if (superfamilies == null) superfamilies = buildSuperfamiliesList();
		return superfamilies;
	}

	protected HashMap<String, ScopDescription> buildSuperfamiliesList() {
		ScopDatabase scop = ScopFactory.getSCOP();
		superfamilies = new HashMap<String, ScopDescription>();
		for (ScopDomain domain : getDomains()) {
			final ScopDescription sf = scop.getScopDescriptionBySunid(domain.getSuperfamilyId());
			superfamilies.put(domain.getScopId(), sf);
		}
		return superfamilies;
	}

}