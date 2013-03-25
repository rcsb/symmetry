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
 * Created on 2013-03-22
 *
 */
package org.biojava3.structure.align.symm.census2;

import java.lang.reflect.Constructor;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

/**
 * A factory for {@link Significance} objects. Most importantly, standardizes definitions of significance across multiple classes. The methods important for this are:
 * <ul>
 * <li>{@link #getForCensus()}, which can be used to determine very liberal significance for running the {@link Census}</li>
 * <li>{@link #getGenerallySymmetric()}, which can be used to determine "general symmetry": both rotational and translational symmetry</li>
 * <li>{@link #getRotationallySymmetric()}, which can be used to determine only rotational symmetry by examining order and angle</li>
 * </ul>
 * These objects, particularly the latter two, may change frequently.
 * Also helps running the {@link Census} with arbitrary Significance objects via {@link CLI}, using the methods:
 * <ul>
 * <li>{@link #SignificanceFactory#fromClass(String)}</li>
 * <li>{@link #SignificanceFactory#fromClass(String, Object...)}</li>
 * <li>{@link #SignificanceFactory#fromMethod(String, String, Object...)}</li>
 * </ul>
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
		return fromMethod(className, methodName, new Object[]{});
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
			throw new IllegalArgumentException("Could not get a Significance object from the method " + methodName + " of "
					+ className, e);
		}
	}

	public static Significance getByConservativeTmScore() {
		return getByTmScore(0.4);
	}

	public static Significance getByLiberalTmScore() {
		return getByTmScore(0.2);
	}

	public static Significance getByMediumTmScore() {
		return getByTmScore(0.3);
	}
	
	public static Significance getByTmScore(final Double cutoff) { // cannot use a primitive double
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

	public static Significance getForCensus() {
		return getByTmScore(0.2);
	}

	public static Significance getGenerallySymmetric() {
		return getByTmScore(0.35);
	}

	public static Significance getRotationallySymmetric() {
		return new Significance() {
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return afpChain.getTMScore() >= 0.3;
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				Axis axis;
				try {
					axis = new Axis(new RotationAxis(afpChain));
				} catch (Exception e) {
					throw new IllegalArgumentException("Could not create a rotation axis for "
							+ protodomain.getString(), e);
				}
				double epsilon = axis.evaluateEpsilon(order);
				if (epsilon > Math.PI / 16) return false;
				return afpChain.getTMScore() >= 0.3 && order >= 2;
			}

			@Override
			public boolean isSignificant(Result result) {
				NumberFormat nf = new DecimalFormat();
				nf.setMaximumFractionDigits(5);
				if (result.getOrder() == null) throw new IllegalArgumentException("Order cannot be null");
				if (result.getAlignment() == null) throw new IllegalArgumentException("Alignment cannot be null");
				final int order = result.getOrder();
				if (order < 2) return false;
				double epsilon = result.getAxis().evaluateEpsilon(order);
				if (epsilon > Math.PI / 16) return false;
				return result.getAlignment().getzScore() >= 3.5;
			}
		};
	}

	public static Significance getUltraLiberal() {
		return getByTmScore(0.0);
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

	private static Class<?>[] params(Object[] args) {
		Class<?>[] classes = new Class<?>[args.length];
		int i = 0;
		for (Object o : args) {
			classes[i] = o.getClass();
			i++;
		}
		return classes;
	}

}
