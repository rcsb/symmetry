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

import org.biojava3.structure.align.symm.census3.CensusResult;


/**
 * A property about a {@link CensusResult}. Usually from {@link CensusResult#getScoreList()}.
 * @author dmyersturnbull
 */
public abstract class CensusResultProperty {

	public abstract double getProperty(CensusResult result) throws PropertyUndefinedException;
	
	public abstract String getName();
	
	@Override
	public String toString() {
		return getName();
	}

	public static CensusResultProperty identity() {
		return new CensusResultProperty() {
			@Override
			public double getProperty(CensusResult result) throws PropertyUndefinedException {
				try {
					return result.getScoreList().getIdentity();
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
	
	public static CensusResultProperty similarity() {
		return new CensusResultProperty() {
			@Override
			public double getProperty(CensusResult result) throws PropertyUndefinedException {
				try {
					return result.getScoreList().getSimilarity();
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

	public static CensusResultProperty order() {
		return new CensusResultProperty() {
			@Override
			public double getProperty(CensusResult result) throws PropertyUndefinedException {
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

	public static CensusResultProperty hasGuessedOrder() {
		return new CensusResultProperty() {
			@Override
			public double getProperty(CensusResult result) throws PropertyUndefinedException {
				if (result.getAxis() == null) return 0;
				return result.getAxis().guessOrder() > 1? 1 : 0;
			}
			@Override
			public String getName() {
				return "hasOrder(angle)";
			}
		};
	}

	public static CensusResultProperty hasOrder() {
		return new CensusResultProperty() {
			@Override
			public double getProperty(CensusResult result) throws PropertyUndefinedException {
				if (result.getGroup() == null || result.getGroup().isAsymmetric()) return 0;
				return 1;
			}
			@Override
			public String getName() {
				return "hasOrder";
			}
		};
	}
	
	public static CensusResultProperty tmScore() {
		return new CensusResultProperty() {
			@Override
			public double getProperty(CensusResult result) throws PropertyUndefinedException {
				try {
					return result.getScoreList().getTmScore();
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

	public static CensusResultProperty epsilon() {
		return new CensusResultProperty() {
			@Override
			public double getProperty(CensusResult result) throws PropertyUndefinedException {
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
