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
 * Created on 2013-02-22
 *
 */
package org.biojava.nbio.structure.align.symm.benchmark.comparison;

import java.util.Random;

import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.symm.census3.CensusResult;


/**
 * A metric that scores symmetry predictions. Has a single method {@link #get(Result)} that determines the quality of a {@link Result}.
 * @author dmyerstu
 *
 */
public abstract class Criterion {

	public abstract double get(CensusResult result) throws NoncomputableCriterionException;

	public abstract String getName();
	
	@Override
	public String toString() {
		return getName();
	}

	public Criterion exp(final double radix) {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				return Math.pow(radix, Criterion.this.get(result));
			}
			@Override
			public String getName() {
				return radix + "^" + Criterion.this.getName();
			}
		};
	}

	public Criterion log() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				return Math.log(Criterion.this.get(result));
			}
			@Override
			public String getName() {
				return "log" + Criterion.this.getName();
			}
		};
	}
	
	public Criterion noFail(final float penalty) {
		return new Criterion() {
			@Override
			public double get(CensusResult result) {
				try {
					return Criterion.this.get(result);
				} catch (NoncomputableCriterionException e) {
					return -penalty;
				}
			}
			@Override
			public String getName() {
				return Criterion.this.getName();
			}
		};
	}
	
	public static Criterion combine(final Criterion a, final Criterion b, final double coeffA, final double coeffB) {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				return coeffA * a.get(result) + coeffB * b.get(result);
			}
			@Override
			public String getName() {
				return coeffA + "*" + a.getName() + " + " + coeffB + "*" + b.getName();
			}
		};
	}

	public static Criterion combineNoFail(final Criterion a, final Criterion b, final double coeffA, final double coeffB) {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				double aa = 0, ba = 0;
				try {
					aa = a.get(result);
				} catch (NoncomputableCriterionException e) {}
				try {
					ba = b.get(result);
				} catch (NoncomputableCriterionException e) {}
				return coeffA * aa + coeffB * ba;
			}
			@Override
			public String getName() {
				return coeffA + "*" + a.getName() + " + " + coeffB + "*" + b.getName();
			}
		};
	}

	public Criterion inverse() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				return -Criterion.this.get(result);
			}
			@Override
			public String getName() {
				return "-" + Criterion.this.getName();
			}
		};
	}

	public static Criterion hasOrder(final float penalty) {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getOrder() < 2) return -penalty;
				return 0;
			}

			@Override
			public String getName() {
				return "hasorder(" + penalty + ")";
			}
		};
	}

	public static Criterion hasOrderLiberal(final float penalty) {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getOrder() < 2) {
					if (result.getAxis() == null) return -penalty;
					if (result.getAxis().guessOrder() == 1) {
						return -penalty;
					}
				}
				return 0;
			}

			@Override
			public String getName() {
				return "hasorderliberal(" + penalty + ")";
			}
		};
	}

	public static Criterion hasOrderByAngle(final float penalty, final double threshold, final int maxOrder) {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				int order = result.getAxis().guessOrder(threshold, maxOrder);
				if (order > 1) return 0;
				return -penalty;
			}

			@Override
			public String getName() {
				return "hasorderbyangle(" + penalty + ")";
			}
		};
	}
	public static Criterion order() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				return result.getOrder();
			}

			@Override
			public String getName() {
				return "order";
			}
		};
	}
	public static Criterion zScore() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getScoreList() == null || result.getScoreList().getzScore() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getScoreList().getzScore();
			}

			@Override
			public String getName() {
				return "Z-score";
			}
		};
	}
	public static Criterion tmScore() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getScoreList() == null || result.getScoreList().getTmScore() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getScoreList().getTmScore();
			}

			@Override
			public String getName() {
				return "TM-score";
			}

		};
	}
	public static Criterion rmsd() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getScoreList() == null || result.getScoreList().getRmsd() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getScoreList().getRmsd();
			}

			@Override
			public String getName() {
				return "RMSD";
			}

		};
	}
	public static Criterion identity() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getScoreList() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getScoreList().getIdentity();
			}

			@Override
			public String getName() {
				return "identity";
			}

		};
	}
	public static Criterion similarity() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getScoreList() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getScoreList().getSimilarity();
			}

			@Override
			public String getName() {
				return "similarity";
			}

		};
	}
	public static Criterion theta() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				return result.getAxis().getAngle();
			}

			@Override
			public String getName() {
				return "theta";
			}

		};
	}
	
	public static Criterion epsilon() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				if (result.getOrder() < 2) throw new NoncomputableCriterionException("The case has a getOrder() of less than 2");
				Double epsilon = result.getAxis().evaluateEpsilon(result.getOrder());
				if (epsilon == null) throw new NoncomputableCriterionException("Could not determine epsilon");
				return (float) ((double) 1.0 / epsilon);
			}

			@Override
			public String getName() {
				return "1/epsilon";
			}

		};
	}
	
	private static volatile Random random = new Random();
	public static Criterion random() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				return random.nextFloat();
			}

			@Override
			public String getName() {
				return "random";
			}

		};
	}
	public static Criterion screw() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getAxis() == null || result.getAxis().getParallel() == null) {
					throw new NoncomputableCriterionException("The case has a null getAxis()");
				}
				return result.getAxis().getParallel();
			}

			@Override
			public String getName() {
				return "screw component";
			}

		};
	}
	
	public static Criterion alignLength() {
		return new Criterion() {
			@Override
			public double get(CensusResult result) throws NoncomputableCriterionException {
				if (result.getScoreList() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getScoreList().getAlignLength();
			}

			@Override
			public String getName() {
				return "aligned length";
			}

		};
	}

}
