package org.biojava3.structure.align.symm.protodomain;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * A symmetry subunit for some SCOP domain. Optionally contains a Structure as {@link #structure}. Responsible for finding the symmetry subunit of a Structure given an appropriately formatted string or an
 * AFPChain and Atom array, and for determining a unique name for this Protodomain.
 * This class is useful for insisting that a protodomain created from an alignment remains within the confines of its domain.
 * To get the protodomain structure string using CE-Symm, do: <code>
 * 		AFPChain afpChain = ceSymm.align(ca1, ca2);
		int order = CeSymm.getSymmetryOrder(afpChain)
		String s = Protodomain.fromSymmetryAlignment(afpChain, ca1, order, atomCache).toString();
 * </code> To get the protdomain structure string by using CE to align a structure against a protodomain of known symmetry (and symmetry order), do: <code>
 * 		AFPChain afpChain = ceSymm.align(ca1, ca2);
		String s = Protodomain.fromReferral(afpChain, ca2, order, atomCache).toString(); // note that we pass in ca2
 * </code>
 * 
 * @author dmyersturnbull
 */
public class Protodomain {

	/**
	 * A helper class to create ResidueRanges without making a mistake.
	 * 
	 * @author dmyersturnbull
	 */
	private static class RangeBuilder {
		private Character currentChain;
		private ResidueNumber currentStart;
		private List<ResidueRange> list = new ArrayList<ResidueRange>();

		void addChain(char chainId) {
			if (currentChain != null)
				throw new IllegalStateException();
			currentChain = chainId;
		}

		void addEnd(ResidueNumber end, int length) {
			if (currentStart == null)
				throw new IllegalStateException();
			list.add(new ResidueRange(currentChain, currentStart, end, length));
			currentStart = null;
			currentChain = null;
		}

		void addStart(ResidueNumber start) {
			if (currentChain == null)
				throw new IllegalStateException();
			if (currentStart != null)
				throw new IllegalStateException();
			currentStart = start;
		}

		List<ResidueRange> getList() {
			return list;
		}
	}

	/**
	 * Client code should avoid calling this method if either fromSymmetryAlignment or fromReferral is more appropriate.
	 * 
	 * @param acaa
	 *            The AFPChain and Atom arrays of the alignment
	 * @param keepStructure
	 *            If set to true, the structure will be recreated.
	 * @param cut
	 *            What the whole protodomain will be "cut" by. This should ordinarily be set to the group order. If it is set to 1, it the whole protodomain will be returned.
	 * @param numBlocks
	 *            The number of blocks in the alignment to consider. Should either be 1 or 2.
	 * @param chainIndex
	 *            The chain in {@code acaa} to create the protodomain from.
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromAfpChain(AFPChain afpChain, Atom[] ca, boolean keepStructure, int cut, int numBlocks,
			int chainIndex, AtomCache cache) throws ProtodomainCreationException {
		ScopDatabase scopInst = ScopFactory.getSCOP();
		ScopDomain scopDomain = scopInst.getDomainByScopID(afpChain.getName2());
		String scopId = scopDomain.getScopId();
		String pdbId = scopDomain.getPdbId();
		// int numAtomsInBlock1Alignment = afpChain.getOptLen()[0];
		// if (numAtomsInBlock1Alignment == 0) throw new ProtodomainCreationException("unknown", scopId);

		// this is a list of insertion-code-independent positions of the groups (where they are in the PDB file)
		Map<ResidueNumber, Integer> posns;
		try {
			Atom[] allAtoms = StructureTools.getAllAtomArray(cache.getStructure(pdbId));
			posns = groupPositions(allAtoms);
		} catch (IOException e1) {
			e1.printStackTrace();
			throw new ProtodomainCreationException("unknown", scopId, e1,
					"Could not get a list of amino acid residue number positions.");
		} catch (StructureException e1) {
			e1.printStackTrace();
			throw new ProtodomainCreationException("unknown", scopId, e1,
					"Could not get a list of amino acid residue number positions.");
		}
		NavigableMap<ResidueNumber, Integer> navMap = orderByAtoms(posns); // we'll need for calculating length

		// since CE-Symm has two blocks, we're going to have to do this in 2 parts:
		// from the start of block 1 to end end of block 1
		// PLUS from the start of block 2 to the end of block 2

		List<String> domainRanges = scopDomain.getRanges();

		int numCaTotal = 0;

		List<ResidueRange> totalRanges = new ArrayList<ResidueRange>();
		// first, we do the first block
		for (int block = 0; block < numBlocks; block++) {
			int numCaInBlock = afpChain.getOptLen()[block];
			numCaTotal += numCaInBlock;
			int[] alignedPositions = afpChain.getOptAln()[block][chainIndex];
			int blockStart = alignedPositions[0];
			int blockEnd = alignedPositions[numCaInBlock - 1];
			Group startGroup = ca[blockStart].getGroup();
			Group endGroup = ca[blockEnd].getGroup();
			ResidueNumber protodomainStartR = startGroup.getResidueNumber();
			ResidueNumber protodomainEndR = endGroup.getResidueNumber();
			int protodomainStart = posns.get(protodomainStartR);
			int protodomainEnd = posns.get(protodomainEndR);
			List<ResidueRange> blockRanges = getBlockRanges(protodomainStart, protodomainEnd, domainRanges,
					protodomainStartR, protodomainEndR, startGroup, endGroup, posns, navMap);
			if (blockRanges == null)
				throw new ProtodomainCreationException("unknown", scopId, "Didn't find any block ranges.");
			totalRanges.addAll(blockRanges);
		}

		List<ResidueRange> cutRanges = cut(totalRanges, posns, numCaTotal / cut); // truncation is okay here

		Protodomain protodomain = new Protodomain(pdbId, cutRanges, cache);

		int length = ResidueRange.calcLength(cutRanges);

		protodomain.length = length; // is the length of the AFPChain
		if (keepStructure)
			try {
				protodomain.createStructure();
			} catch (IOException e) {
				e.printStackTrace();
				throw new ProtodomainCreationException(protodomain.toString(), scopId, e,
						"Could not create the structure, which was required.");
			} catch (StructureException e) {
				e.printStackTrace();
				throw new ProtodomainCreationException(protodomain.toString(), scopId, e,
						"Could not create the structure, which was required.");
			}

		return protodomain;

	}

	/**
	 * Creates a protodomain from an alignment between a Protodomain of known symmetry and some other structure.
	 * @param afpChain
	 *            The <em>referred</em> structure (result) must be in position 2, and the <em>referring</em> structure (query) in position 1.
	 * @param ca
	 *            An array of Atoms of the result
	 * @param cut
	 *            What the whole protodomain will be "cut" by. This should ordinarily be set to the group order of the query. If it is set to 1, it the whole protodomain will be returned.
	 * @param cache
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromReferral(AFPChain afpChain, Atom[] ca, int cut, AtomCache cache)
			throws ProtodomainCreationException {
		try {
			return fromAfpChain(afpChain, ca, false, cut, 1, 1, cache);
		} catch (Exception e) { // IOException | StructureException | ArrayIndexOutOfBoundsException |
			// NullPointerException
			e.printStackTrace();
			throw new ProtodomainCreationException("unknown", afpChain.getName2(), e);
		}
	}

	/**
	 * @param string
	 *            A string of the form: pdbId.chain_start-end,chain_start-end, ..., chain_start-end.
	 * @param scopDomain
	 * @param keepStructure
	 *            If set to true, the structure will be recreated.
	 * @param cut
	 *            What the whole protodomain will be "cut" by. This should ordinarily be set to the group order. If it is set to 1, it the whole protodomain will be returned.
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromString(String string, ScopDomain scopDomain, boolean keepStructure, int cut,
			AtomCache cache) throws ProtodomainCreationException {

		String pdbId = scopDomain.getPdbId();

		List<ResidueRange> list = ResidueRange.parseMultiple(string.substring(string.indexOf('.') + 1));

		Protodomain protodomain = new Protodomain(pdbId, list, cache);
		if (keepStructure) {
			try {
				protodomain.createStructure();
			} catch (IOException e) {
				e.printStackTrace();
				throw new ProtodomainCreationException(protodomain.toString(), pdbId, e);
			} catch (StructureException e) {
				e.printStackTrace();
				throw new ProtodomainCreationException(protodomain.toString(), pdbId, e);
			}
		}

		return protodomain;

	}

	/**
	 * @param string
	 *            A string of the form: pdbId.chain_start-end,chain_start-end, ..., chain_start-end.
	 * @param scopId
	 * @param keepStructure
	 *            If set to true, the structure will be recreated.
	 * @param cut
	 *            What the whole protodomain will be "cut" by. This should ordinarily be set to the group order. If it is set to 1, it the whole protodomain will be returned.
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromString(String string, String scopId, boolean keepStructure, int cut, AtomCache cache)
			throws ProtodomainCreationException, IllegalArgumentException {
		return fromString(string, ScopFactory.getSCOP().getDomainByScopID(scopId), keepStructure, cut, cache);
	}

	/**
	 * Creates a Protodomain from an AFPChain returned by CE-Symm.
	 * @param afpChain
	 * @param ca
	 *            An array of Atoms of the result
	 * @param cut
	 *            What the whole protodomain will be "cut" by. This should ordinarily be set to the group order. If it is set to 1, it the whole protodomain will be returned.
	 * @param cache
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromSymmetryAlignment(AFPChain afpChain, Atom[] ca, int cut, AtomCache cache)
			throws ProtodomainCreationException {
		return fromAfpChain(afpChain, ca, false, cut, 2, 1, cache);
	}

	/**
	 * A way to work around insertion codes.
	 * 
	 * @param ca
	 *            An array of Amino acid C-alpha Atoms from a single Structure
	 * @return A HashMap in which each ResidueNumber of an amino acid Group found in {@code ca} is mapped to its exact position in the ATOM fields of the PDB file
	 */
	public static Map<ResidueNumber, Integer> groupPositions(Atom[] atoms) {
		HashMap<ResidueNumber, Integer> hm = new HashMap<ResidueNumber, Integer>(); // TODO remove me! change back to hashmap
		for (int i = 0; i < atoms.length; i++) {
			Group g = atoms[i].getGroup();
			ResidueNumber rn = g.getResidueNumber();
			if (g.getType().equals(GroupType.AMINOACID)) {
				if (!hm.containsKey(rn)) {
					hm.put(rn, i + 1);
				}
			}
		}
		return hm;
	}

	/**
	 * From a list of residue ranges {@code rrs} and a number of residues wanted {@code numResiduesWant} (<em>which includes any alignment gaps</em>), returns a new list of residue ranges that
	 * contains the first {@code numResiduesWant} residues from {@code rrs}.
	 * 
	 * @param rrs
	 * @param groupPositions
	 *            A map returned from groupPositions()
	 * @param numResiduesWant
	 * @return
	 */
	private static List<ResidueRange> cut(List<ResidueRange> rrs, Map<ResidueNumber, Integer> groupPositions,
			int numResiduesWant) {

		List<ResidueRange> part = new ArrayList<ResidueRange>(); // the parts we want to keep

		// we'll need these in a bit
		NavigableMap<ResidueNumber, Integer> navMap = orderByAtoms(groupPositions);

		int numResiduesHave = 0;

		outer: for (ResidueRange rr : rrs) {

			// note that getLength() here DOES INCLUDE gaps (and it should)
			if (numResiduesHave + rr.getLength() <= numResiduesWant) {
				// no problem: we can fit the entire residue range in
				part.add(rr);

			} else if (numResiduesHave <= numResiduesWant) {
				// okay, so we can only fit some of rr in (and then we'll end)
				// we'll insert rr from rr.getStart() to rr.getStart() + numResiduesWant - numResiduesHave
				// unfortunately, this is harder because of insertion codes
				for (Map.Entry<ResidueNumber, Integer> entry : navMap.entrySet()) { // this is in
					// insertion-code-independent order
					if (entry.getKey().equals(rr.getStart())) {
						ResidueNumber endResidueNumber = entry.getKey(); // where we want to end this last residue range
						// okay, we got to the end of rr
						// now we need to go back by numResiduesWant - numResiduesTraversed
						for (int j = 0; j < numResiduesWant - numResiduesHave; j++) {
							endResidueNumber = navMap.higherEntry(endResidueNumber).getKey();
						}
						part.add(new ResidueRange(rr.getChain(), rr.getStart(), endResidueNumber, numResiduesWant
								- numResiduesHave));
						break outer; // we filled up numResiduesTraversed, so we won't be able to add further residue ranges
					}
				}
			}
			numResiduesHave += rr.getLength();
		}
		return part;
	}

	/**
	 * Finds the list of residue ranges that corresponds to the protodomain between {@code protodomainStartR} and {@code protodomainEndR} <em>within the SCOP domain identified by {@code domainRanges}.
	 * 
	 * @param protodomainStart
	 * @param protodomainEnd
	 * @param domainRanges
	 * @param protodomainStartR
	 * @param protodomainEndR
	 * @param blockStartGroup
	 * @param blockEndGroup
	 * @param posns
	 * @param navMap
	 * @return
	 */
	private static List<ResidueRange> getBlockRanges(int protodomainStart, int protodomainEnd,
			List<String> domainRanges, ResidueNumber protodomainStartR, ResidueNumber protodomainEndR,
			Group blockStartGroup, Group blockEndGroup, Map<ResidueNumber, Integer> posns,
			NavigableMap<ResidueNumber, Integer> navMap) {

		/*
		 * Okay, things get annoying here, since we need to handle:
		 * 1. domain insertions, AND
		 * 2. heteroatom Groups inside this domain that are past the last amino acid
		 * If we simply take the protodomain's start and its end (in this block), we could include:
		 * 1. domains inserted between, AND
		 * 2. residues in other domains between the last amino acid in the correct domain and the last heteroatom (especially a water molecule) in the correct domain
		 * These are bad. So: we use SCOP's list of ranges for our domain, making sure that our protodomain does not extend outside
		 * its boundaries.
		 * To handle #2, we -could- get the last amino acid (for the end residue; the first aa for the start), but this is problematic for modified amino acids. We don't do this.
		 * Some heteroatoms also have alpha carbons.
		 * But that's not all. What about multi-chain domains? We need to handle those.
		 * Then things can get even MORE cumbersome thanks to insertion codes. d1qdmc1
		 * has an insertion code in one of its SCOP ranges: C:1S-104S (Residue numbers in ATOMs go: 247 1S 2S ... 103S 104S 248)
		 */

		if (protodomainStart > protodomainEnd)
			return null; // it's okay if our protodomain doesn't exist in this block

		RangeBuilder rangeBuilder = new RangeBuilder();
		for (String domainRange : domainRanges) {
			// a range is of the form: A:05-117 (chain:start-end) OR just A:
			ResidueNumber domainStartR, domainEndR;
			int domainStart, domainEnd;
			if (domainRange.length() > 1) {
				// chain:start-end
				String[] myParts = domainRange.substring(2).split("-");
				domainStartR = ResidueNumber.fromString(myParts[0]);
				domainEndR = ResidueNumber.fromString(myParts[1]);
			} else {
				// no actual numbers; just the chain
				domainStartR = blockStartGroup.getResidueNumber();
				domainEndR = blockEndGroup.getResidueNumber();
			}
			String chainId = domainRange.substring(0, 1);
			domainStartR.setChainId(chainId);
			domainEndR.setChainId(chainId);

			// these are the insertion-code-independent positions
			domainStart = posns.get(domainStartR);
			domainEnd = posns.get(domainEndR);

			if (domainStart >= domainEnd)
				return null; // in this case, there's something wrong with the SCOP definition

			char chain = domainRange.charAt(0);
			rangeBuilder.addChain(chain);

			if (domainEnd <= protodomainStart) { // the domain part starts and ends before the start of the protodomain
				continue; // well, obviously we're not using that part of the domain
			} else if (domainStart <= protodomainStart) { // the domain part starts before the protodomain starts but
				// ends after the protodomain starts
				rangeBuilder.addStart(protodomainStartR); // protodomain start
				if (domainEnd <= protodomainEnd) { // the domain part ends before the protodomain ends
					rangeBuilder.addEnd(domainEndR, residuesBetween(domainEnd, protodomainStart, navMap, chain)); // domain end
				} else { // the domain part ends after the protodomain ends
					rangeBuilder.addEnd(protodomainEndR,
							residuesBetween(protodomainEnd, protodomainStart, navMap, chain)); // protodomain end
				}
			} else if (domainStart > protodomainEnd) { // domain part ends before the protodomain starts
				continue; // if we knew the domain parts are ordered, we could break
			} else if (domainStart > protodomainStart) { // the domain part starts after the protodomain starts but ends
				// after it starts
				rangeBuilder.addStart(domainStartR); // domain start
				if (domainEnd <= protodomainEnd) {
					rangeBuilder.addEnd(domainEndR, residuesBetween(domainEnd, domainStart, navMap, chain)); // domain end
				} else {
					rangeBuilder.addEnd(protodomainEndR, residuesBetween(protodomainEnd, domainStart, navMap, chain)); // protodomain end
				}
			} else { // this can't happen
				throw new RuntimeException(); // might as well catch a bug immediately
			}
		}

		List<ResidueRange> list = rangeBuilder.getList();
		return list;
	}

	/**
	 * @param original
	 *            A map of amino acid Group ATOM record positions, in any order
	 * @return A map of amino acid Group ATOM record positions, ordered by value (position).
	 */
	private static NavigableMap<ResidueNumber, Integer> orderByAtoms(Map<ResidueNumber, Integer> original) {
		Comparator<ResidueNumber> vc = new ValueComparator<ResidueNumber, Integer>(original);
		NavigableMap<ResidueNumber, Integer> navMap = new TreeMap<ResidueNumber, Integer>(vc);
		navMap.putAll(original);
		// for (Map.Entry<ResidueNumber, Integer> entry : navMap.entrySet()) System.out.println(entry.getKey() + "\t\t\t" + entry.getValue());
		return navMap;
	}

	/**
	 * @param positionA
	 * @param positionB
	 * @param orderedPosns
	 * @param chain
	 * @return The number of amino acids between {@code positionA} and {@code positionB} in the ATOM records.
	 */
	private static int residuesBetween(int positionA, int positionB, NavigableMap<ResidueNumber, Integer> orderedPosns,
			char chain) {
		int positionStart, positionEnd;
		if (positionA <= positionB) {
			positionStart = positionA;
			positionEnd = positionB;
		} else {
			positionStart = positionB;
			positionEnd = positionA;
		}
		int count = 0;
		for (Map.Entry<ResidueNumber, Integer> entry : orderedPosns.entrySet()) {
			if (entry.getKey().getChainId().charAt(0) == chain) {
				if (entry.getValue() == positionStart) {
					count = 0;
				}
				if (entry.getValue() == positionEnd)
					return count;
				count++;
			}
		}
		return -1;
	}

	private AtomCache cache;

	private Integer length;

	private final List<ResidueRange> list;

	private final String pdbId;

	private final String string;

	private Structure structure;

	public Protodomain(String pdbId, List<ResidueRange> list, AtomCache cache) {
		this.pdbId = pdbId;
		this.list = list;
		this.cache = cache;
		StringBuilder sb = new StringBuilder();
		Iterator<ResidueRange> iter = list.iterator();
		while (iter.hasNext()) {
			ResidueRange t = iter.next();
			sb.append(t);
			if (iter.hasNext())
				sb.append(",");
		}
		string = sb.toString();
	}

	/**
	 * Builds the structure of this Protodomain.
	 * 
	 * @throws IOException
	 * @throws StructureException
	 * @see #getStructure()
	 */
	public void createStructure() throws IOException, StructureException {
		if (structure == null) {
			structure = cache.getStructure(toString());
			structure.setName(toString());
		}
	}

	public AtomCache getCache() {
		return cache;
	}

	/**
	 * @return The list of residue ranges in this Protodomain.
	 */
	public List<ResidueRange> getList() {
		return list;
	}

	public String getPdbId() {
		return pdbId;
	}

	/**
	 * <strong>This method has not been tested yet! Do not use.</strong>
	 * 
	 * @return The number of amino acids in this Protodomain, including gaps, or null if the length is unknown. The length will be unknown if and only if this Protodomain was created using a string
	 *         rather than an AFPChain. AFPChain.
	 * 
	 */
	public Integer getRangeLength() {
		return length;
	}

	/**
	 * Same as {@link #toString()}.
	 * 
	 * @return
	 */
	public String getString() {
		return toString();
	}

	/**
	 * @return The Structure of this Protodomain, or null if it has not been created. To create the structure, call {@link #createStructure()}.
	 */
	public Structure getStructure() {
		return structure;
	}

	/**
	 * Nullifies the structure of this Protodomain (presumably to save memory).
	 * 
	 * @see #getStructure()
	 */
	public void removeStructure() {
		structure = null;
	}

	public void setCache(AtomCache cache) {
		this.cache = cache;
	}

	@Override
	public String toString() {
		return pdbId + "." + string;
	}

}
