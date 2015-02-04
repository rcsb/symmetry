package org.biojava.nbio.structure.align.symm.census3.stats;

import org.biojava.nbio.structure.align.symm.census3.CensusResult;

/**
 * 
 * @author dmyersturnbull
 *
 */
public class CensusResultPropertyFactory {

	private static abstract class Callback<T> {
		abstract T getProperty(CensusResult result);
	}
	
	private static <T extends Number> CensusResultProperty<T> build(final Callback<T> callback, final String name) {
		return new CensusResultProperty<T>() {
			@Override
			public String getName() {
				return name;
			}
			@Override
			public T getProperty(CensusResult result) throws PropertyUndefinedException {
				try {
					return callback.getProperty(result);
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String toString() {
				return name;
			}
		};
	}
	
	public static CensusResultProperty<Double> epsilon() {
		return build(new Callback<Double>() {
			@Override
			Double getProperty(CensusResult result) {
				return result.getAxis().evaluateEpsilon(result.getOrder());
			}
		}, "epsilon");
	}

	public static CensusResultProperty<Byte> hasAnyOrder() {
		return build(new Callback<Byte>() {
			@Override
			Byte getProperty(CensusResult result) {
				if (result.getGroup() != null) {
					if (!result.getGroup().isAsymmetric()) return (byte) 1.0;
				}
				if (result.getAxis() != null) {
					if (result.getAxis().guessOrder() > 1) return (byte) 1.0;
				}
				return (byte) 0.0;
			}
		}, "[order>1]");
	}

	public static CensusResultProperty<Float> identity() {
		return build(new Callback<Float>() {
			@Override
			Float getProperty(CensusResult result) {
				return result.getScoreList().getIdentity();
			}
		}, "identity");
	}
	
	public static CensusResultProperty<Float> similarity() {
		return build(new Callback<Float>() {
			@Override
			Float getProperty(CensusResult result) {
				return result.getScoreList().getSimilarity();
			}
		}, "similarity");
	}

	public static CensusResultProperty<Float> tmScore() {
		return build(new Callback<Float>() {
			@Override
			Float getProperty(CensusResult result) {
				return result.getScoreList().getTmScore();
			}
		}, "TM-score");
	}
	
}
