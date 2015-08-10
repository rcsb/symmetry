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
package org.biojava.nbio.structure.align.symm.census3;

import java.lang.reflect.Constructor;

/**
 * A factory for {@link CensusSignificance} objects. Most importantly, standardizes definitions of significance across
 * multiple classes. The methods important for this are:
 * <ul>
 * <li>{@link #forCensus()}, which can be used to determine very liberal significance for running the {@link Census}</li>
 * <li>{@link #generallySymmetric()}, which can be used to determine "general symmetry": both rotational and
 * translational symmetry</li>
 * <li>{@link #rotationallySymmetric()}, which can be used to determine only rotational symmetry by examining order and
 * angle</li>
 * </ul>
 * These objects, particularly the latter two, may change frequently. Also helps running the {@link Census} with
 * arbitrary Significance objects via {@link CensusCLI}, using the methods:
 * <ul>
 * <li>{@link #SignificanceFactory#fromClass(String)}</li>
 * <li>{@link #SignificanceFactory#fromClass(String, Object...)}</li>
 * <li>{@link #SignificanceFactory#fromMethod(String, String, Object...)}</li>
 * </ul>
 * 
 * @author dmyerstu
 */
public class CensusSignificanceFactory {

	public static CensusSignificance and(final CensusSignificance a, final CensusSignificance b) {
		return new CensusSignificance() {
			@Override
			public boolean isSignificant(CensusResult result) {
				return a.isSignificant(result) && b.isSignificant(result);
			}
		};
	}

	public static CensusSignificance or(final CensusSignificance a, final CensusSignificance b) {
		return new CensusSignificance() {
			@Override
			public boolean isSignificant(CensusResult result) {
				return a.isSignificant(result) || b.isSignificant(result);
			}
		};
	}

	public static CensusSignificance fromClass(String className) {
		try {
			return (CensusSignificance) Class.forName(className).newInstance();
		} catch (Exception e) {
			throw new IllegalArgumentException("Could not instantiate " + className, e);
		}
	}

	public static CensusSignificance fromClass(String className, Object... args) {
		try {
			Class<?> c = Class.forName(className);
			Constructor<?> con = c.getConstructor(params(args));
			return (CensusSignificance) con.newInstance(args);
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
	public static CensusSignificance fromMethod(String className, String methodName) {
		return fromMethod(className, methodName, new Object[] {});
	}

	/**
	 * Gets a Significance object returned from the method {@code methodName} in the class {@code className}.
	 * 
	 * @param className
	 *            A <em>fully-qualified</em> class name; e.g.
	 *            org.biojava.nbio.structure.align.symm.census2.SignificanceFactory. If null, defaults to this this class's
	 *            name.
	 * @param methodName
	 * @param args
	 * @return
	 */
	public static CensusSignificance fromMethod(String className, String methodName, Object... args) {
		if (className == null) className = CensusSignificanceFactory.class.getName();
		try {
			return (CensusSignificance) Class.forName(className).getMethod(methodName, params(args)).invoke(null, args);
		} catch (Exception e) {
			throw new IllegalArgumentException("Could not get a Significance object from the method " + methodName
					+ " of " + className, e);
		}
	}

	public static CensusSignificance not(final CensusSignificance a) {
		return new CensusSignificance() {
			@Override
			public boolean isSignificant(CensusResult result) {
				return !a.isSignificant(result);
			}
		};
	}

	public static CensusSignificance notForCeSymmOrd() {
		return not(and(tmScore(0.4), order()));
	}

	public static CensusSignificance forCeSymmOrd() {
		return and(tmScore(0.4), order());
	}

	public static CensusSignificance forCeSymmTm() {
		return tmScore(0.4);
	}

	public static CensusSignificance forPublishedSymD10() {
		return symdZScore(10);
	}

	public static CensusSignificance forPublishedSymD8() {
		return symdZScore(8);
	}

	public static CensusSignificance forUnpublishedSymD() {
		return tmScore(0.47);
	}

	public static CensusSignificance order() {
		return new CensusSignificance() {
			@Override
			public boolean isSignificant(CensusResult result) {
				if (result.getOrder() > 1) return true;
				if (result.getAxis() != null) {
					return result.getAxis().guessOrder(1.0 * Math.PI / 180, 8) > 1;
				}
				return false;
			}
		};
	}

	public static CensusSignificance tmScore(final Double cutoff) { // cannot use a primitive double
		return new CensusSignificance() {
			@Override
			public boolean isSignificant(CensusResult result) {
				return result.getScoreList().getTmScore() >= cutoff;
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

	private static CensusSignificance symdZScore(final double cutoff) {
		return new CensusSignificance() {
			@Override
			public boolean isSignificant(CensusResult result) {
				return result.getScoreList().getzScore() >= cutoff;
			}
		};
	}

	public static CensusSignificance ultraLiberal() {
		return new CensusSignificance() {
			@Override
			public boolean isSignificant(CensusResult result) {
				return true;
			}
		};
	}

}
