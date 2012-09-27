package org.biojava3.structure.align.symm.protodomain;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.ResidueNumber;

/**
 * A chain, a start residue, and an end residue. May also store a length value; because of insertion codes, this length is not necessarily end-start.
 * @author dmyersturnbull
 *
 */
public class ResidueRange {
	private final char chain;
	private final ResidueNumber start;
	private final ResidueNumber end;
	private final Integer length;
	
	public String toString() {
		return chain + "_" + start + "-" + end;
	}
	
	/**
	 * @return The number of residues in this ResidueRange, including alignment gaps. This value will be null if and only if this ResidueRange was created with a null length.
	 */
	public Integer getLength() {
		return length;
	}
	public ResidueRange(char chain, ResidueNumber start, ResidueNumber end, Integer length) {
		this.chain = chain;
		this.start = start;
		this.end = end;
		this.length = length;
	}
	public char getChain() {
		return chain;
	}
	public ResidueNumber getStart() {
		return start;
	}
	public ResidueNumber getEnd() {
		return end;
	}
	/**
	 * Calculates the combined number of residues of the ResidueRanges in {@code rrs}, <em>given that each ResidueRange has a length calculated</em>.
	 * The value, if calculated, <em>will include any alignment gaps</em>.
	 * @param rrs A list of ResidueRanges
	 * @return The combined length
	 * @throws IllegalArgumentException If a ResidueRange's length is null
	 * @see #getLength()
	 */
	public static int calcLength(List<ResidueRange> rrs) {
		int l = 0;
		for (ResidueRange rr : rrs) {
			if (rr.getLength() == null) throw new IllegalArgumentException("A ResidueRange does not have a length.");
			l += rr.getLength();
		}
		return l;
	}
	
	/**
	 * @param s A string of the form chain_start-end
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static ResidueRange parse(String s) {
		char chain = s.charAt(0);
		String[] parts = s.substring(2).split("-");
		ResidueNumber start = ResidueNumber.fromString(parts[0]);
		ResidueNumber end = ResidueNumber.fromString(parts[1]);
		return new ResidueRange(chain, start, end, null);
	}
	
	/**
	 * @param s A string of the form chain_start-end,chain_start-end, ...
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static List<ResidueRange> parseMultiple(String s) {
		String[] parts = s.split(",");
		List<ResidueRange> list = new ArrayList<ResidueRange>(parts.length);
		for (String part : parts) {
			list.add(parse(part));
		}
		return list;
	}
}
