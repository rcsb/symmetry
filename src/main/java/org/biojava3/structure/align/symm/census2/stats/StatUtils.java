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

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Static utilities for statistics.
 * @author dmyerstu
 */
public class StatUtils {

	public static final String NEWLINE;
	private static final int MAX_FRACTION_DIGITS = 3;

	private static NumberFormat nf = new DecimalFormat();

	static {
		NEWLINE = System.getProperty("line.separator");
	}
	static {
		nf.setMaximumFractionDigits(MAX_FRACTION_DIGITS);
		nf.setMinimumFractionDigits(1);
	}

	public static String formatD(double d) {
		return nf.format(d);
	}

	public static String formatP(double p) {
		return nf.format(p*100.0) + "%";
	}
	
	public static <T> void plus(Map<T, Double> map, T key, double value) {
		if (!map.containsKey(key)) map.put(key, 0.0);
		map.put(key, map.get(key) + value);
	}

	public static <T> void plus(Map<T, Integer> map, T key) {
		plus(map, key, 1);
	}

	public static <T> void plus(Map<T, Integer> map, T key, int value) {
		if (!map.containsKey(key)) map.put(key, 0);
		map.put(key, map.get(key) + value);
	}

	public static <T,V> void plusSet(Map<T, Set<V>> map, T key, V value) {
		if (!map.containsKey(key)) map.put(key, new HashSet<V>());
		map.get(key).add(value);
	}

	public static <T> void plusD(Map<T, Double> map, T key) {
		plus(map, key, 1);
	}

	public static <T> void plusD(Map<T, Double> map, T key, double value) {
		if (!map.containsKey(key)) map.put(key, 0.0);
		map.put(key, map.get(key) + value);
	}

}
