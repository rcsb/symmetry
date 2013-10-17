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
package org.biojava3.structure.align.symm.census2.stats;

import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava3.structure.align.symm.census2.Result;

/**
 * A unit by which to report statistics. Particularly: family, superfamily, or fold.
 * @author dmyerstu
 */
public abstract class Grouping {

	@Override
	public abstract String toString();
	
	public abstract String group(Result result);
	
	public abstract String group(ScopDomain domain);

	public static Grouping byName(String name) {
		try {
			return (Grouping) Grouping.class.getMethod(name, new Class[]{}).invoke(null, new Object[]{});
		} catch (Exception e) {
			throw new RuntimeException("Could not get a Grouping object from the method " + name, e);
		}
	}
	
	public static Grouping domain() {
		return new Grouping() {
			@Override
			public String group(Result result) {
				return result.getScopId();
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
	public static Grouping superfamily() {
		return new Grouping() {
			@Override
			public String group(Result result) {
				String[] parts = result.getClassification().split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + result.getScopId());
				return parts[0] + "." + parts[1] + "." + parts[2];
			}

			@Override
			public String toString() {
				return "superfamily";
			}
			@Override
			public String group(ScopDomain domain) {
				String[] parts = domain.getClassificationId().split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + domain.getScopId());
				return parts[0] + "." + parts[1] + "." + parts[2];
			}
		};
	}

	public static Grouping family() {
		return new Grouping() {
			@Override
			public String group(Result result) {
				String[] parts = result.getClassification().split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + result.getScopId());
				return parts[0] + "." + parts[1] + "." + parts[2] + "." + parts[3];
			}

			@Override
			public String toString() {
				return "family";
			}
			@Override
			public String group(ScopDomain domain) {
				String[] parts = domain.getClassificationId().split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + domain.getScopId());
				return parts[0] + "." + parts[1] + "." + parts[2] + "." + parts[3];
			}
		};
	}

	public static Grouping fold() {
		return new Grouping() {
			@Override
			public String group(Result result) {
				String[] parts = result.getClassification().split("\\.");
				if (parts.length < 2) throw new IllegalArgumentException("Classification id is invalid for " + result.getScopId());
				return parts[0] + "." + parts[1];
			}

			@Override
			public String toString() {
				return "fold";
			}
			@Override
			public String group(ScopDomain domain) {
				String[] parts = domain.getClassificationId().split("\\.");
				if (parts.length < 3) throw new IllegalArgumentException("Classification id is invalid for " + domain.getScopId());
				return parts[0] + "." + parts[1];
			}
		};
	}

}
