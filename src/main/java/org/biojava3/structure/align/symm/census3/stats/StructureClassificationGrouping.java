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
 * Created on 2013-03-10
 *
 */
package org.biojava3.structure.align.symm.census3.stats;

import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census3.CensusResult;

/**
 * A unit by which to report statistics. Particularly: family, superfamily, or fold.
 * @author dmyersturnbull
 */
public abstract class StructureClassificationGrouping {

	@Override
	public abstract String toString();
	
	public abstract String group(CensusResult result);
	
	public abstract String group(ScopDomain domain);

	public static StructureClassificationGrouping byName(String name) {
		try {
			return (StructureClassificationGrouping) StructureClassificationGrouping.class.getMethod(name, new Class[]{}).invoke(null, new Object[]{});
		} catch (Exception e) {
			throw new RuntimeException("Could not get a Grouping object from the method " + name, e);
		}
	}
	
	public static StructureClassificationGrouping domain() {
		return new StructureClassificationGrouping() {
			@Override
			public String group(CensusResult result) {
				return result.getId();
			}
			@Override
			public String toString() {
				return "domain";
			}
			@Override
			public String group(ScopDomain domain) {
				return domain.getScopId();
			}
		};
	}
	public static StructureClassificationGrouping superfamily() {
		return new StructureClassificationGrouping() {
			@Override
			public String group(CensusResult result) {
				String classification = ScopFactory.getSCOP().getDomainByScopID(result.getId()).getClassificationId();
				String[] parts = classification.split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + result.getId());
				return parts[0] + "." + parts[1] + "." + parts[2];
			}

			@Override
			public String toString() {
				return "superfamily";
			}
			@Override
			public String group(ScopDomain domain) {
				String classification = ScopFactory.getSCOP().getDomainByScopID(domain.getScopId()).getClassificationId();
				String[] parts = classification.split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + domain.getScopId());
				return parts[0] + "." + parts[1] + "." + parts[2];
			}
		};
	}

	public static StructureClassificationGrouping family() {
		return new StructureClassificationGrouping() {
			@Override
			public String group(CensusResult result) {
				String classification = ScopFactory.getSCOP().getDomainByScopID(result.getId()).getClassificationId();
				String[] parts = classification.split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + result.getId());
				return parts[0] + "." + parts[1] + "." + parts[2] + "." + parts[3];
			}

			@Override
			public String toString() {
				return "family";
			}
			@Override
			public String group(ScopDomain domain) {
				String classification = ScopFactory.getSCOP().getDomainByScopID(domain.getScopId()).getClassificationId();
				String[] parts = classification.split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + domain.getScopId());
				return parts[0] + "." + parts[1] + "." + parts[2] + "." + parts[3];
			}
		};
	}

	public static StructureClassificationGrouping fold() {
		return new StructureClassificationGrouping() {
			@Override
			public String group(CensusResult result) {
				String classification = ScopFactory.getSCOP().getDomainByScopID(result.getId()).getClassificationId();
				String[] parts = classification.split("\\.");
				if (parts.length < 2) throw new IllegalArgumentException("Classification id is invalid for " + result.getId());
				return parts[0] + "." + parts[1];
			}

			@Override
			public String toString() {
				return "fold";
			}
			@Override
			public String group(ScopDomain domain) {
				String[] parts = domain.getClassificationId().split("\\.");
				if (parts.length < 2) throw new IllegalArgumentException("Classification id is invalid for " + domain.getScopId());
				return parts[0] + "." + parts[1];
			}
		};
	}

	public static StructureClassificationGrouping clas() {
		return new StructureClassificationGrouping() {
			@Override
			public String group(CensusResult result) {
				String classification = ScopFactory.getSCOP().getDomainByScopID(result.getId()).getClassificationId();
				String[] parts = classification.split("\\.");
				if (parts.length < 1) throw new IllegalArgumentException("Classification id is invalid for " + result.getId());
				return parts[0];
			}

			@Override
			public String toString() {
				return "class";
			}
			@Override
			public String group(ScopDomain domain) {
				String[] parts = domain.getClassificationId().split("\\.");
				if (parts.length < 1) throw new IllegalArgumentException("Classification id is invalid for " + domain.getScopId());
				return parts[0];
			}
		};
	}

}
