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

import org.biojava3.structure.align.symm.census2.Result;

/**
 * 
 * @author dmyerstu
 */
public abstract class Property {

	public abstract double getProperty(Result result) throws PropertyUndefinedException;
	
	public abstract String getName();
	
	@Override
	public String toString() {
		return getName();
	}

	public static Property identity() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getAlignment().getIdentity();
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "identity";
			}
		};
	}
	
	public static Property similarity() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getAlignment().getSimilarity();
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "similarity";
			}
		};
	}

	public static Property order() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getOrder();
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "order";
			}
		};
	}

	public static Property hasOrder() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				if (result.getOrder() == null) return 0;
				return 1;
			}
			@Override
			public String getName() {
				return "order";
			}
		};
	}
	
	public static Property tmScore() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getAlignment().getTmScore();
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "TM-score";
			}
		};
	}

	public static Property epsilon() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getAxis().evaluateEpsilon(result.getOrder());
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "epsilon";
			}
		};
	}
	
}
