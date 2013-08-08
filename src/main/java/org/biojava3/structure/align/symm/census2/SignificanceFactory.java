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
 * Created on 2013-03-22
 */
package org.biojava3.structure.align.symm.census2;

import java.lang.reflect.Constructor;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

/**
 * A factory for {@link Significance} objects. Most importantly, standardizes definitions of significance across
 * multiple classes. The methods important for this are:
 * <ul>
 * <li>{@link #forCensus()}, which can be used to determine very liberal significance for running the {@link Census}</li>
 * <li>{@link #generallySymmetric()}, which can be used to determine "general symmetry": both rotational and
 * translational symmetry</li>
 * <li>{@link #rotationallySymmetric()}, which can be used to determine only rotational symmetry by examining order and
 * angle</li>
 * </ul>
 * These objects, particularly the latter two, may change frequently. Also helps running the {@link Census} with
 * arbitrary Significance objects via {@link CLI}, using the methods:
 * <ul>
 * <li>{@link #SignificanceFactory#fromClass(String)}</li>
 * <li>{@link #SignificanceFactory#fromClass(String, Object...)}</li>
 * <li>{@link #SignificanceFactory#fromMethod(String, String, Object...)}</li>
 * </ul>
 * 
 * @author dmyerstu
 */
public class SignificanceFactory {

	public static Significance and(final Significance a, final Significance b) {
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return a.isPossiblySignificant(afpChain) && b.isPossiblySignificant(afpChain);
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				return a.isSignificant(protodomain, order, angle, afpChain)
						&& b.isSignificant(protodomain, order, angle, afpChain);
			}

			@Override
			public boolean isSignificant(Result result) {
				return a.isSignificant(result) && b.isSignificant(result);
			}
		};
	}

	public static Significance asymmetric() {
		return asymmetric(0.4);
	}

	public static Significance asymmetric(final double tmScore) {
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return afpChain.getTMScore() < tmScore;
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				if (order < 2) return true;
				return afpChain.getTMScore() < tmScore;
			}

			@Override
			public boolean isSignificant(Result result) {
				if (result.getOrder() == null || result.getOrder() < 2) return true;
				return result.getAlignment().getTmScore() < tmScore;
			}
		};
	}

	public static Significance bin1() {
		return binned(0.4, 0.5);
	}

	public static Significance bin2() {
		return binned(0.5, 0.6);
	}

	public static Significance bin3() {
		return binned(0.6, 1);
	}

	public static Significance binned(final double start, final double end) {
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return afpChain.getTMScore() >= start && afpChain.getTMScore() < end;
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				if (order < 2) return false;
				return afpChain.getTMScore() >= start && afpChain.getTMScore() < end;
			}

			@Override
			public boolean isSignificant(Result result) {
				if (result.getOrder() == null || result.getOrder() < 2) return false;
				return result.getAlignment().getTmScore() >= start && result.getAlignment().getTmScore() < end;
			}
		};
	}

	public static Significance conservative() {
		return rotationallySymmetric(0.5, Math.PI);
	}

	public static Significance forCensus() {
		return SignificanceFactory.generallySymmetric();
	}

	public static Significance forCeSymmOrd() {
		return rotationallySymmetricSmart();
	}

	public static Significance forCeSymmTm() {
		return rotationallySymmetricSmart();
	}

	public static Significance forPublishedSymD10() {
		return zScore(10);
	}

	public static Significance forPublishedSymD8() {
		return zScore(8);
	}

	public static Significance forUnpublishedSymD() {
		return tmScore(0.47);
	}

	public static Significance fromClass(String className) {
		try {
			return (Significance) Class.forName(className).newInstance();
		} catch (Exception e) {
			throw new IllegalArgumentException("Could not instantiate " + className, e);
		}
	}

	public static Significance fromClass(String className, Object... args) {
		try {
			Class<?> c = Class.forName(className);
			Constructor<?> con = c.getConstructor(params(args));
			return (Significance) con.newInstance(args);
		} catch (Exception e) {
			throw new IllegalArgumentException("Could not instantiate " + className, e);
		}
	}

/**
		 * Calls {@link #fromMethod(String, String, Object...) with no args for the method to be called.
		 * @param className
		 * @param methodName
		 * @return
		 */
	public static Significance fromMethod(String className, String methodName) {
		return fromMethod(className, methodName, new Object[] {});
	}

	/**
	 * Gets a Significance object returned from the method {@code methodName} in the class {@code className}.
	 * 
	 * @param className
	 *            A <em>fully-qualified</em> class name; e.g.
	 *            org.biojava3.structure.align.symm.census2.SignificanceFactory. If null, defaults to this this class's
	 *            name.
	 * @param methodName
	 * @param args
	 * @return
	 */
	public static Significance fromMethod(String className, String methodName, Object... args) {
		if (className == null) className = SignificanceFactory.class.getName();
		try {
			return (Significance) Class.forName(className).getMethod(methodName, params(args)).invoke(null, args);
		} catch (Exception e) {
			throw new IllegalArgumentException("Could not get a Significance object from the method " + methodName
					+ " of " + className, e);
		}
	}

	public static Significance generallySymmetric() {
		return tmScore(0.4);
	}

	public static Significance liberal() {
		return rotationallySymmetric(0.3, Double.MAX_VALUE);
	}

	public static Significance liberalTmScore() {
		return tmScore(0.3);
	}

	public static Significance mediumTmScore() {
		return tmScore(0.4);
	}

	public static Significance not(final Significance a) {
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return !a.isPossiblySignificant(afpChain);
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				return !a.isSignificant(protodomain, order, angle, afpChain);
			}

			@Override
			public boolean isSignificant(Result result) {
				return !a.isSignificant(result);
			}
		};
	}

	public static Significance notRotationallySymmetricSmart() {
		return SignificanceFactory.not(or(rotationallySymmetric(), rotationallySymmetricWithAngle()));
	}

	public static Significance or(final Significance a, final Significance b) {
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return a.isPossiblySignificant(afpChain) || b.isPossiblySignificant(afpChain);
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				return a.isSignificant(protodomain, order, angle, afpChain)
						|| b.isSignificant(protodomain, order, angle, afpChain);
			}

			@Override
			public boolean isSignificant(Result result) {
				return a.isSignificant(result) || b.isSignificant(result);
			}
		};
	}

	public static Significance rotationallySymmetric() {
		return rotationallySymmetric(0.4, Double.MAX_VALUE);
	}

	public static Significance rotationallySymmetric(final double threshold, final double epsilonThreshold) {

		return new Significance() {

			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return afpChain.getTMScore() >= threshold;
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				Axis axis;
				try {
					axis = new Axis(new RotationAxis(afpChain));
				} catch (Exception e) {
					return false;
				}
				double epsilon = axis.evaluateEpsilon(order);
				if (epsilon > epsilonThreshold) return false;
				return afpChain.getTMScore() >= threshold && order >= 2;
			}

			@Override
			public boolean isSignificant(Result result) {
				NumberFormat nf = new DecimalFormat();
				nf.setMaximumFractionDigits(5);
				if (result.getOrder() == null) return false;
				if (result.getAlignment() == null) return false;
				final int order = result.getOrder();
				if (order < 2) return false;
				double epsilon = result.getAxis().evaluateEpsilon(order);
				if (epsilon > epsilonThreshold) return false;
				return result.getAlignment().getTmScore() >= threshold;
			}
		};
	}

	public static Significance rotationallySymmetricSmart() {
		return or(rotationallySymmetric(), rotationallySymmetricWithAngle());
	}

	public static Significance rotationallySymmetricWithAngle() {
		return rotationallySymmetricWithAngle(0.4, 8, 1.0 * Math.PI / 180);
	}

	public static Significance rotationallySymmetricWithAngle(final double threshold, final int maxOrder,
			final double deviationThreshold) {
		return new Significance() {

			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return afpChain.getTMScore() >= threshold;
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				Axis axis;
				try {
					axis = new Axis(new RotationAxis(afpChain));
				} catch (Exception e) {
					return false;
				}
				int theOrder = axis.guessOrder(deviationThreshold, maxOrder);
				return afpChain.getTMScore() >= threshold && theOrder > 1;
			}

			@Override
			public boolean isSignificant(Result result) {
				if (result.getAxis() == null) return false;
				if (result.getAlignment() == null) return false;
				int theOrder = result.getAxis().guessOrder(deviationThreshold, maxOrder);
				return result.getAlignment().getTmScore() >= threshold && theOrder > 1;
			}
		};
	}

	public static Significance superConservative() {
		return rotationallySymmetric(0.7, Math.PI);
	}

	public static Significance tmScore(final Double cutoff) { // cannot use a primitive double
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return afpChain.getTMScore() >= cutoff;
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				return afpChain.getTMScore() >= cutoff;
			}

			@Override
			public boolean isSignificant(Result result) {
				return result.getAlignment().getTmScore() >= cutoff;
			}
		};
	}

	public static Significance ultraLiberal() {
		return tmScore(0.0);
	}

	public static Significance veryConservative() {
		return rotationallySymmetric(0.6, Math.PI);
	}

	private static Class<?>[] params(Object[] args) {
		Class<?>[] classes = new Class<?>[args.length];
		int i = 0;
		for (Object o : args) {
			classes[i] = o.getClass();
			i++;
		}
		return classes;
	}

	private static Significance zScore(final double cutoff) {
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return afpChain.getProbability() >= cutoff;
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				return afpChain.getProbability() >= cutoff;
			}

			@Override
			public boolean isSignificant(Result result) {
				return result.getAlignment().getzScore() >= cutoff;
			}
		};
	}

}
