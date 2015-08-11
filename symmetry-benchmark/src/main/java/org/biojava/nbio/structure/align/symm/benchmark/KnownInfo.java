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
 * Created on 2013-03-03
 */
package org.biojava.nbio.structure.align.symm.benchmark;

import java.io.Serializable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A known (true) space group and corresponding order of symmetry. Can differentiate between different types of
 * symmetry.
 * 
 * @author dmyerstu
 */
public class KnownInfo implements Serializable, Comparable<KnownInfo> {
	private static final long serialVersionUID = -2667699023747790086L;
	private String group;

	public static int getOrderFromGroup(String string) {
		int order = 1;
		Pattern pattern = Pattern.compile("([\\d])$");
		Matcher matcher = pattern.matcher(string);
		boolean found = matcher.find();
		if (found && matcher.groupCount() == 1) {
			order = Integer.parseInt(matcher.group(1));
		}
		return order;
	}

	public KnownInfo() {

	}

	public KnownInfo(String group) {
		super();
		this.group = group;
	}

	@Override
	public int compareTo(KnownInfo o) {
		if (equals(o)) return 0;
		if (getOrder() < o.getOrder()) return -1;
		if (getOrder() > o.getOrder()) return 1;
		return group.compareTo(o.getGroup());
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		KnownInfo other = (KnownInfo) obj;
		if (group == null) {
			if (other.group != null) return false;
		} else if (!group.equals(other.group)) return false;
		return true;
	}

	public String getGroup() {
		return group;
	}

	public int getOrder() {
		return getOrderFromGroup(group);
	}

	/**
	 * @return Whether this info has <em>cyclic but not dihedral</em> symmetry.
	 */
	public boolean hasCyclicSymmetry() {
		return !isAsymmetric() && group.contains("C");
	}

	/**
	 * @return Whether this info has dihedral symmetry.
	 */
	public boolean hasDihedralSymmetry() {
		return !isAsymmetric() && group.contains("D");
	}

	/**
	 * @return Whether this info has <em>rotational or true helical</em> symmetry and its order of symmetry is an even
	 *         number.
	 */
	public boolean hasEvenOrderSymmetry() {
		return (hasRotationalSymmetry() || hasTrueHelicalSymmetry()) && getOrder() % 2 == 0;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (group == null ? 0 : group.hashCode());
		return result;
	}

	/**
	 * @return Whether this info <em>has helical but not</em> rotational symmetry.
	 *         <em>This includes superhelical and non-integral order helical symmetry.</em>
	 */
	public boolean hasHelicalSymmetry() {
		return !isAsymmetric() && group.contains("H");
	}

	/**
	 * @return Whether this info has helical symmetry for which the symmetry operation must be applied a fractional
	 *         number of times to result in the trivial alignment.
	 *         <em>This includes superhelical and non-integral order helical symmetry.</em>
	 */
	public boolean hasNonIntegralOrderSymmetry() {
		return group.equals("NIH");
	}

	/**
	 * @return Whether this info has <em>rotational or true helical</em> symmetry and its order of symmetry is an odd
	 *         number.
	 */
	public boolean hasOddOrderSymmetry() {
		return (hasRotationalSymmetry() || hasTrueHelicalSymmetry()) && getOrder() % 2 == 1;
	}

	/**
	 * @return Whether this info has true rotational symmetry.
	 */
	public boolean hasRotationalSymmetry() {
		return !isAsymmetric() && group.contains("C") || group.contains("D");
	}

	/**
	 * @return Whether this info has a "bent" or curved helical symmetry, such as that of a screw that is bent around in
	 *         a curve.
	 */
	public boolean hasSuperhelicalSymmetry() {
		return group.startsWith("SH");
	}

	/**
	 * @return Whether this info <em>has translational but not</em> rotational symmetry.
	 */
	public boolean hasTranslationalSymmetry() {
		return !isAsymmetric() && group.contains("R");
	}

	/**
	 * @return Whether this info <em>has <strong>true</strong> helical symmetry</em>. This <em>does not</em> include
	 *         superhelical or non-integral order helical symmetry.
	 */
	public boolean hasTrueHelicalSymmetry() {
		if (group.startsWith("H")) {
			System.err.println(group);
		}
		return group.startsWith("H");
	}

	/**
	 * @return Whether this info has no symmetry of any kind.
	 */
	public boolean isAsymmetric() {
		return group.equals("C1");
	}

	public void setGroup(String group) {
		this.group = group;
	}

	@Override
	public String toString() {
		return group;
	}

}