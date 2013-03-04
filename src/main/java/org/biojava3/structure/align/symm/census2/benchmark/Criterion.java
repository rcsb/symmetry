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
package org.biojava3.structure.align.symm.census2.benchmark;

import java.util.Random;

import org.biojava3.structure.align.symm.census2.Result;

/**
 * A classifier for the quality of CE-Symm results. Has a method {@link #get(Case)} that determines the quality of a {@link Result}. Can also {@link #hasSymmetry(Result) decide} whether a result is symmetric.
 * @author dmyerstu
 *
 * @param <T> The numerical type of the quality
 */
public abstract class Criterion<T extends Number> {

	public static Criterion<Float> inverseF(final Criterion<Float> criterion) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				return -criterion.get(result);
			}

			@Override
			public String getName() {
				return "-(" + criterion.getName() + ")";
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return !criterion.hasSymmetry(result);
			}
		};
	}

	public static Criterion<Integer> inverseI(final Criterion<Integer> criterion) {
		return new Criterion<Integer>() {
			@Override
			public Integer get(Result result) throws NoncomputableCriterionException {
				return -criterion.get(result);
			}

			@Override
			public String getName() {
				return "-(" + criterion.getName() + ")";
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return !criterion.hasSymmetry(result);
			}
		};
	}

	public static Criterion<Float> combineFF(final Criterion<Float> a, final Criterion<Float> b, final float coeffA, final float coeffB) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				return coeffA * a.get(result) + coeffB * b.get(result);
			}

			@Override
			public String getName() {
				return coeffA + "*(" + a.getName() + ") + " + coeffB + " *(" + b.getName() + ")";
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return a.hasSymmetry(result) && b.hasSymmetry(result);
			}
		};
	}

	public static Criterion<Float> combineFI(final Criterion<Float> a, final Criterion<Integer> b, final float coeffA, final float coeffB) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				return coeffA * a.get(result) + coeffB * b.get(result);
			}

			@Override
			public String getName() {
				return coeffA + "*(" + a.getName() + ") + " + coeffB + " *(" + b.getName() + ")";
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return a.hasSymmetry(result) && b.hasSymmetry(result);
			}
		};
	}

	public static Criterion<Integer> combineII(final Criterion<Integer> a, final Criterion<Integer> b, final float coeffA, final float coeffB) {
		return new Criterion<Integer>() {
			@Override
			public Integer get(Result result) throws NoncomputableCriterionException {
				return (int) (coeffA * a.get(result) + coeffB * b.get(result));
			}

			@Override
			public String getName() {
				return coeffA + "*(" + a.getName() + ") + " + coeffB + " *(" + b.getName() + ")";
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return a.hasSymmetry(result) && b.hasSymmetry(result);
			}
		};
	}
	
	public abstract T get(Result result) throws NoncomputableCriterionException;
	
	public abstract String getName();
	
	public abstract boolean hasSymmetry(Result result) throws NoncomputableCriterionException;

	public static Criterion<Float> zScore(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAlignment() == null || result.getAlignment().getzScore() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getAlignment().getzScore();
			}

			@Override
			public String getName() {
				return "Z-score>" + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Float> tmScore(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAlignment() == null || result.getAlignment().getTmScore() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getAlignment().getTmScore();
			}

			@Override
			public String getName() {
				return "TM-score>" + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Float> rmsd(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAlignment() == null || result.getAlignment().getRmsd() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getAlignment().getRmsd();
			}

			@Override
			public String getName() {
				return "RMSD>" + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Float> identity(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getAlignment().getIdentity();
			}

			@Override
			public String getName() {
				return "identity>" + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Float> similarity(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getAlignment().getSimilarity();
			}

			@Override
			public String getName() {
				return "similarity>" + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Float> theta(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				return result.getAxis().getTheta();
			}

			@Override
			public String getName() {
				return "theta>" + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Float> thetaIsCorrect(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				if (result.getOrder() == null) throw new NoncomputableCriterionException("The case has a null getOrder()");
				float theta = result.getAxis().getTheta();
				float shouldBe = (float) (2*Math.PI / result.getOrder());
				return Math.abs(shouldBe - theta);
			}

			@Override
			public String getName() {
				return "2pi/o - theta > " + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Float> epsilon(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				if (result.getOrder() == null) throw new NoncomputableCriterionException("The case has a null getOrder()");
				if (result.getOrder() < 2) throw new NoncomputableCriterionException("The case has a getOrder() of less than 2");
				Double epsilon = result.getAxis().evaluateEpsilon(result.getOrder());
				if (epsilon == null) throw new NoncomputableCriterionException("Could not determine epsilon");
				return (float) ((double) 1.0 / epsilon);
			}

			@Override
			public String getName() {
				return "1/epsilon > " + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	private static volatile Random random = new Random();
	public static Criterion<Float> random(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				return random.nextFloat();
			}

			@Override
			public String getName() {
				return "random (" + (threshold*100.0f) + "%)";
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Float> screw(final float threshold) {
		return new Criterion<Float>() {
			@Override
			public Float get(Result result) throws NoncomputableCriterionException {
				if (result.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				System.out.println(result.getAxis().getScrew());
				return result.getAxis().getScrew();
			}

			@Override
			public String getName() {
				return "screw component > " + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	public static Criterion<Integer> alignLength(final float threshold) {
		return new Criterion<Integer>() {
			@Override
			public Integer get(Result result) throws NoncomputableCriterionException {
				if (result.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return result.getAlignment().getAlignLength();
			}

			@Override
			public String getName() {
				return "aligned length > " + threshold;
			}

			@Override
			public boolean hasSymmetry(Result result) throws NoncomputableCriterionException {
				return get(result) > threshold;
			}
		};
	}
	
}
