package org.biojava3.structure.align.symm.census2.benchmark;


public abstract class Criterion<T extends Number> {

	public abstract T get(Case c) throws NoncomputableCriterionException;

	public static Criterion<Float> zScore() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return c.getAlignment().getzScore();
			}
		};
	}
	public static Criterion<Float> tmScore() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return c.getAlignment().getTmScore();
			}
		};
	}
	public static Criterion<Float> rmsd() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return c.getAlignment().getRmsd();
			}
		};
	}
	public static Criterion<Float> identity() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return c.getAlignment().getIdentity();
			}
		};
	}
	public static Criterion<Float> similarity() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return c.getAlignment().getSimilarity();
			}
		};
	}
	public static Criterion<Float> theta() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				return c.getAxis().getTheta();
			}
		};
	}
	public static Criterion<Float> thetaIsCorrect() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				if (c.getOrder() == null) throw new NoncomputableCriterionException("The case has a null getOrder()");
				float theta = c.getAxis().getTheta();
				float shouldBe = (float) (2*Math.PI / c.getOrder());
				return Math.abs(shouldBe - theta);
			}
		};
	}
	public static Criterion<Float> epsilon() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				if (c.getOrder() == null) throw new NoncomputableCriterionException("The case has a null getOrder()");
				if (c.getOrder() < 2) throw new NoncomputableCriterionException("The case has a getOrder() of less than 2");
				Double epsilon = c.getAxis().evaluateEpsilon(c.getOrder());
				if (epsilon == null) throw new NoncomputableCriterionException("Could not determine epsilon");
				return (float) ((double) epsilon);
			}
		};
	}
	public static Criterion<Float> screw() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) throws NoncomputableCriterionException {
				if (c.getAxis() == null) throw new NoncomputableCriterionException("The case has a null getAxis()");
				return c.getAxis().getScrew();
			}
		};
	}
	public static Criterion<Integer> alignLength() {
		return new Criterion<Integer>() {
			@Override
			public Integer get(Case c) throws NoncomputableCriterionException {
				if (c.getAlignment() == null) throw new NoncomputableCriterionException("The case has a null getAlignment()");
				return c.getAlignment().getAlignLength();
			}
		};
	}
	
}
