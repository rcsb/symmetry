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
 * Created on 2013-03-03
 *
 */
package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.Serializable;

/**
 * A known (true) group and order of symmetry.
 * @author dmyerstu
 */
public class KnownInfo implements Serializable, Comparable<KnownInfo> {
	private static final long serialVersionUID = -2667699023747790086L;
	private String group;
	private int order;
	public String getGroup() {
		return group;
	}
	public int getOrder() {
		return order;
	}
	public KnownInfo(String group, int order) {
		super();
		this.group = group;
		this.order = order;
	}
	public KnownInfo() {
		
	}
	
	public void setGroup(String group) {
		this.group = group;
	}
	public void setOrder(int order) {
		this.order = order;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((group == null) ? 0 : group.hashCode());
		return result;
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
	
	@Override
	public String toString() {
		return "KnownInfo [group=" + group + ", order=" + order + "]";
	}
	@Override
	public int compareTo(KnownInfo o) {
		if (equals(o)) return 0;
		if (order < o.getOrder()) return -1;
		if (order > o.getOrder()) return 1;
		return group.compareTo(o.getGroup());
	}
}