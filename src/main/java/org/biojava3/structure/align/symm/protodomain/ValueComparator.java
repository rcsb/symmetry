package org.biojava3.structure.align.symm.protodomain;

import java.util.Comparator;
import java.util.Map;

/**
 * A map that is sorted by its values.
 * @author dmyersturnbull
 *
 * @param <T> The key type
 * @param <V> The value type
 */
public class ValueComparator<T, V extends Comparable<V>> implements Comparator<T> {

	private Map<T, V> map;

	public ValueComparator(Map<T, V> map) {
		this.map = map;
	}

	@Override
	public int compare(T o1, T o2) {
		return map.get(o1).compareTo(map.get(o2));
	}

}
