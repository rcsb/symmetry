package org.biojava3.structure.align.symm.benchmark.comparison;

/**
 * Something that can decide whether two space groups are equivalent, and whether two orders of rotational symmetry are
 * equivalent.
 * 
 * @author dmyerstu
 */
public interface GroupGuesser {
	boolean hasEquivalentGroup(String a, String b);

	boolean hasEquivalentOrder(int a, int b);
}
