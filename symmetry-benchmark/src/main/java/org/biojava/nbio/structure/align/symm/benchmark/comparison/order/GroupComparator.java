package org.biojava.nbio.structure.align.symm.benchmark.comparison.order;

/**
 * Something that can decide whether two space groups are equivalent, and whether two orders of rotational symmetry are
 * equivalent.
 * 
 * @author dmyerstu
 */
public interface GroupComparator {
	boolean hasEquivalentGroup(String a, String b);

	boolean hasEquivalentOrder(int a, int b);
}
