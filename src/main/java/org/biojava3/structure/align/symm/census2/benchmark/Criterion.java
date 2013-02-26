package org.biojava3.structure.align.symm.census2.benchmark;

public abstract class Criterion<T extends Number> {

	public abstract T get(Case c);

	public static Criterion<Float> zScore() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) {
				return c.getAlignment().getzScore();
			}
		};
	}
	public static Criterion<Float> tmScore() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) {
				return c.getAlignment().getTmScore();
			}
		};
	}
	public static Criterion<Float> rmsd() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) {
				return c.getAlignment().getRmsd();
			}
		};
	}
	public static Criterion<Float> identity() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) {
				return c.getAlignment().getIdentity();
			}
		};
	}
	public static Criterion<Float> similarity() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) {
				return c.getAlignment().getSimilarity();
			}
		};
	}
	public static Criterion<Float> theta() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) {
				return c.getAxis().getTheta();
			}
		};
	}
	public static Criterion<Float> thetaIsCorrect() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) {
				float theta = c.getAxis().getTheta();
				float shouldBe = (float) (2*Math.PI / c.getOrder());
				return Math.abs(shouldBe - theta);
			}
		};
	}
	public static Criterion<Float> screw() {
		return new Criterion<Float>() {
			@Override
			public Float get(Case c) {
				return c.getAxis().getScrew();
			}
		};
	}
	public static Criterion<Integer> alignLength() {
		return new Criterion<Integer>() {
			@Override
			public Integer get(Case c) {
				return c.getAlignment().getAlignLength();
			}
		};
	}
	
}
