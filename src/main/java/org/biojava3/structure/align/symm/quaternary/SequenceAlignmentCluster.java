package org.biojava3.structure.align.symm.quaternary;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;

public class SequenceAlignmentCluster {
	private double sequenceIdentityThreshold = 1.0;
	private List<UniqueSequenceList> uniqueSequenceList = new ArrayList<UniqueSequenceList>();
	private List<Atom[]> alignedCAlphaAtoms = null;
	private List<Atom[]> alignedCBetaAtoms = null;
	private boolean modified = true;

	public SequenceAlignmentCluster (double sequenceIdentityThreshold) {
		this.sequenceIdentityThreshold = sequenceIdentityThreshold;
	}
	
	public boolean addChain(Atom[] cAlphaAtoms, String chainId) {
		// check if new chain exactly maches an exisiting sequence, then add it
		if (exactMatch(cAlphaAtoms, chainId)) {
			modified = true;
			return true;
		}
		// if not an exact match, try a sequence similarity match
		if (similarityMatch(cAlphaAtoms, chainId)) {
			modified = true;
			return true;
		}
		return false;
	}

	public void addUniqueSequenceList(UniqueSequenceList sequenceList) {
		uniqueSequenceList.add(sequenceList);
		modified = true;
	}

	public int getSequenceCount() {
		int count = 0;
		for (UniqueSequenceList list: uniqueSequenceList) {
			count += list.getChainCount();
		}
		return count;
	}
	
	public int getSequenceAlignmentLength() {
		return getReferenceResidueIndices().size();
	}
	
	public List<String> getChainIds() {
		List<String> ids = new ArrayList<String>();
		for (UniqueSequenceList list: uniqueSequenceList) {
			for (int i = 0, n = list.getChainCount(); i < n; i++) {
				ids.add(list.getChainId(i));
			}
		}
		return ids;
	}
	
	public List<Atom[]> getAlignedCalphaAtoms() {
		run();
		return alignedCAlphaAtoms;
	}
	
	public List<Atom[]> getAlignedCBetaAtoms() {
		run();
		return alignedCBetaAtoms;
	}

	private void run() {
		if (modified) {
			alignedCBetaAtoms = null;
			alignedCAlphaAtoms = null;
			createAlignedCAlphaAtoms();
			createAlignedCBetaAtoms();	
			modified = false;
		}
	}

	public List<Atom[]> getCbetaAtoms() {
		return new ArrayList<Atom[]>();
	}

	public List<Integer> getSubunitEquivalenceClasses() {
		return new ArrayList<Integer>();
	}

	public String toString() {
		StringBuilder builder = new StringBuilder();
		for (UniqueSequenceList u: uniqueSequenceList) {
			builder.append(u.toString());
			builder.append("\n");
		}
		return builder.toString();
	}
	
	private boolean exactMatch(Atom[] cAlphaAtoms, String chainId) {
		for (UniqueSequenceList u: uniqueSequenceList) {
			if (u.addChain(cAlphaAtoms, chainId)) {
				List<Integer> align = new ArrayList<Integer>(cAlphaAtoms.length);
				// add identity sequence alignment
				for (int i = 0; i < cAlphaAtoms.length; i++) {
					align.add(i);
				}
				System.out.println("SequenceAlignmentCluster: found exact match");
				return true;
			}
		}
		return false;
	}

	private boolean similarityMatch(Atom[] cAlphaAtoms, String chainId) {
		UniqueSequenceList u = uniqueSequenceList.get(0);
		AFPChain afp = alignPair(u.getReferenceChain(), cAlphaAtoms);
		if (afp == null) {
			return false;
		}
		double identity = afp.getIdentity();
		if (identity < sequenceIdentityThreshold) {	
			return false;
		}

		int[][][] alignment = afp.getOptAln();
		if (alignment != null) {
			UniqueSequenceList seqList = new UniqueSequenceList(cAlphaAtoms, chainId);
			List<Integer> align1 = new ArrayList<Integer>(alignment[0][0].length);
			for (Integer a1: alignment[0][0]) {
				align1.add(a1);
			}
			seqList.setAlignment1(align1);

			List<Integer> align2 = new ArrayList<Integer>(alignment[0][1].length);
			for (Integer a2: alignment[0][1]) {
				align2.add(a2);
			}
			seqList.setAlignment2(align2);
			System.out.println("Sequence alignment: ");
			System.out.println(align1);
			System.out.println(align2);
			addUniqueSequenceList(seqList);
			return true;
		}
		return false;
	}

	private AFPChain alignPair(Atom[] ca1Seq, Atom[] ca2Seq) {
//		System.out.println("SequenceAlignmentCluster: aligning sequences");
		SmithWaterman3Daligner aligner = new SmithWaterman3Daligner();
		AFPChain afp = null;
		try {
			afp = aligner.align(ca1Seq, ca2Seq);
		} catch (StructureException e) {
			e.printStackTrace();
			return afp;
		} 
		return afp;
	}
	
	private void createAlignedCAlphaAtoms() {
		List<Integer> indices = getReferenceResidueIndices();
		System.out.println("In common: " + indices);
		alignedCAlphaAtoms = new ArrayList<Atom[]>();
		for (UniqueSequenceList u: uniqueSequenceList) {
			List<Integer> alignment1 = u.getAlignment1();
			List<Integer> alignment2 = u.getAlignment2();
			List<Integer> alignmentIndices = new ArrayList<Integer>();
			for (int i = 0; i < alignment1.size(); i++) {
				int a1 = alignment1.get(i);
				if (indices.contains(a1)) {
					alignmentIndices.add(alignment2.get(i));
				}
			}
			for (int i = 0; i < u.getChainCount(); i++) {
				Atom[] unalignedAtoms = u.getChain(i);
				Atom[] alignedAtoms = new Atom[alignmentIndices.size()];
				System.out.println("Unaligned: " + unalignedAtoms.length);
				System.out.println("Aligned: " + alignedAtoms.length);
				for (int j = 0; j < alignedAtoms.length; j++) {
					alignedAtoms[j] = unalignedAtoms[alignmentIndices.get(j)];
				}
				System.out.println("Calpha: " + alignedAtoms.length);
				alignedCAlphaAtoms.add(alignedAtoms);
			}
		}
	}
	
	private void createAlignedCBetaAtoms() {
		List<Atom[]>cbetas = new ArrayList<Atom[]>(alignedCAlphaAtoms.size());
		for (Atom[] cas: alignedCAlphaAtoms) {
		    Atom[] cbs = new Atom[cas.length];
			for (int i = 0, n = cas.length; i < n; i++) {
				try {
					Atom b = cas[i].getGroup().getAtomByPDBname(" CB ");
	//				System.out.println("Found cb: " + b.getName());
					if (b != null) {
						cbs[i] = b;
					}
				} catch (StructureException e) {;
				}
			}
			cbetas.add(cbs);
		}
		
		// create a set that contains the indices of atoms
		// that are in common among all chains
		Set<Integer> inCommon = new HashSet<Integer>();
		Atom[] cb = cbetas.get(0);
		for (int i = 0; i < cb.length; i++) {
			if (cb[i] != null) {
				inCommon.add(i);
			}
		}
		
		Set<Integer> occurances = new HashSet<Integer>();
		// note j must start at 1 here!
		for (int j = 1; j < cbetas.size(); j++) {
			cb = cbetas.get(j);
			for (int i = 0; i < cb.length; i++) {
				if (cb[i] != null) {
					occurances.add(i);
				}
			}
			inCommon.retainAll(occurances);
			occurances.clear();
		}
		
//		System.out.println("CBs in common: " + inCommon);
		// copy only atoms that are in common in all chains
		alignedCBetaAtoms = new ArrayList<Atom[]>();
		for (int j = 0; j < cbetas.size(); j++) {
			cb = cbetas.get(j);
			Atom[] cbatoms = new Atom[inCommon.size()];
			int index = 0;
			for (int i = 0; i < cb.length; i++) {
			   if (inCommon.contains(i)) {
//				   if (cb[i] == null) {
//					   System.out.println("Adding null: " + index);
//				   }
				   cbatoms[index] = cb[i];
				   index++;
			   }
			}
//			System.out.println("index: " + index);
			alignedCBetaAtoms.add(cbatoms);
			System.out.println("Cbeta: " + cbatoms.length);
		}
	}
	
	private List<Integer> getReferenceResidueIndices() {
		List<Integer> indices = new ArrayList<Integer>(uniqueSequenceList.get(0).getAlignment1());
		for (UniqueSequenceList u: uniqueSequenceList) {
           indices.retainAll(u.getAlignment1());
		}
		return indices;
	}
}
