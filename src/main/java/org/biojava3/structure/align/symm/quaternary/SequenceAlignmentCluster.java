package org.biojava3.structure.align.symm.quaternary;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3DParameters;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

public class SequenceAlignmentCluster {
	private QuatSymmetryParameters parameters = null;
	private List<UniqueSequenceList> uniqueSequenceList = new ArrayList<UniqueSequenceList>();
	private List<Atom[]> alignedCAlphaAtoms = null;
//	private List<Atom[]> alignedCBetaAtoms = null;
	private int alignmentLength = 0;
	private boolean modified = true;

	public SequenceAlignmentCluster (QuatSymmetryParameters parameters) {
		this.parameters = parameters;
	}
	
	/**
	 * Returns true if the passed in sequence matches this sequence alignment cluster within the
	 * within the sequence identity threshold and alignment fraction threshold
	 * @param sequence
	 * @return true if sequence matches cluster
	 */
	public boolean isSequenceMatch(String sequence) {
		if (uniqueSequenceList.size() == 0) {
			return false;
		}
		
		UniqueSequenceList u = uniqueSequenceList.get(0);
		String referenceSequence = u.getSeqResSequence();
		ProteinSequence s1 = new ProteinSequence(referenceSequence);
		if (s1.getLength() == 0) {
			return false;
		}
		
		ProteinSequence s2 = new ProteinSequence(sequence);
		if (s2.getLength() == 0) {
			return false;
		}
		
		if (referenceSequence.equals(sequence)) {
			return true;
		}
		
		SmithWaterman3DParameters params = new SmithWaterman3DParameters();
		GapPenalty penalty = new SimpleGapPenalty();
		penalty.setOpenPenalty(params.getGapOpen());
		penalty.setExtensionPenalty(params.getGapExtend());
		SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();
		PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman =
			Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

		SequencePair<ProteinSequence, AminoAcidCompound> pair = smithWaterman.getPair();
//		System.out.println("Len1: " + referenceSequence.length() + " - " + s1.getLength());
//		System.out.println("Len2: " + sequence.length() + " - " + s2.getLength());
//		System.out.println("Length: " + pair.getLength() + "idenities: " + pair.getNumIdenticals());
		if (pair.getLength() == 0)
			return false;
		double sequenceIdentity = (double)pair.getNumIdenticals()/pair.getLength();
		double alignmentLengthFraction = (double)pair.getLength()/Math.max(s1.getLength(), s2.getLength());
//		System.out.println("Ref seq. identity: " + sequenceIdentity + " fraction: " + alignmentLengthFraction);
		return sequenceIdentity >= parameters.getSequenceIdentityThreshold() && alignmentLengthFraction >= parameters.getAlignmentFractionThreshold();
	}
	
	public boolean addChain(Atom[] cAlphaAtoms, String chainId, int modelNumber, String sequence) {	
		// check if new chain exactly matches an existing sequence, then add it to the cluster
		if (exactMatch(cAlphaAtoms, chainId, modelNumber)) {
			modified = true;
			return true;
		}
		// check is SEQRES are identical
		if (identityMatch(cAlphaAtoms, chainId, modelNumber, sequence)) {
			modified = true;
			return true;
		}
		// if not an exact match, try a sequence similarity match
		if (similarityMatch(cAlphaAtoms, chainId, modelNumber, sequence)) {
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
		run();
		return alignmentLength;
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
	
	public List<Integer> getModelNumbers() {
		List<Integer> numbers = new ArrayList<Integer>();
		for (UniqueSequenceList list: uniqueSequenceList) {
			for (int i = 0, n = list.getChainCount(); i < n; i++) {
				numbers.add(list.getModelNumber(i));
			}
		}
		return numbers;
	}
	
	public List<Atom[]> getAlignedCalphaAtoms() {
		run();
		return alignedCAlphaAtoms;
	}
	
//	public List<Atom[]> getAlignedCBetaAtoms() {
//		run();
//		return alignedCBetaAtoms;
//	}
	
	public String toString() {
		StringBuilder builder = new StringBuilder();
		for (UniqueSequenceList u: uniqueSequenceList) {
			builder.append(u.toString());
			builder.append("\n");
		}
		return builder.toString();
	}

	private void run() {
		if (modified) {
//			alignedCBetaAtoms = null;
			alignedCAlphaAtoms = null;
			createAlignedCAlphaAtoms();
			createAlignedCBetaAtoms();	
			modified = false;
		}
	}
	
	private boolean exactMatch(Atom[] cAlphaAtoms, String chainId, int modelNumber) {
		for (UniqueSequenceList u: uniqueSequenceList) {
			if (u.isMatch(cAlphaAtoms)) {
			    u.addChain(cAlphaAtoms, chainId, modelNumber);
				List<Integer> align = new ArrayList<Integer>(cAlphaAtoms.length);
				// add identity sequence alignment
				for (int i = 0; i < cAlphaAtoms.length; i++) {
					align.add(i);
				}
//				System.out.println("SequenceAlignmentCluster: found exact match");
				return true;
			}
		}
		return false;
	}
	
	private boolean identityMatch(Atom[] cAlphaAtoms, String chainId, int modelNumber, String sequence) {
		UniqueSequenceList u = uniqueSequenceList.get(0);

		// check for 100% identity match of reference sequence
		String refSequence = u.getSeqResSequence();
		boolean seqMatch = refSequence.equals(sequence);

		// if reference (SEQRES) sequences match 100%, 
		// find alignment of atom sequences by Smith-Waterman alignment
		if (seqMatch) {
			List<Integer> alig1 = new ArrayList<Integer>();
			List<Integer> alig2 = new ArrayList<Integer>();
			Atom[] referenceAtoms = u.getReferenceChain();
			int inCommon = alignByAtomSequence(referenceAtoms, cAlphaAtoms, alig1, alig2);
//			System.out.println("in common: "  + inCommon);

			UniqueSequenceList seqList = new UniqueSequenceList(cAlphaAtoms, chainId, modelNumber, sequence);
			seqList.setAlignment1(alig1);
			seqList.setAlignment2(alig2);
			//			System.out.println(alig1);
			//			System.out.println(alig2);
			addUniqueSequenceList(seqList);
			return true;
		}

		return false;
	}

	private boolean similarityMatch(Atom[] cAlphaAtoms, String chainId, int modelNumber, String sequence) {
		UniqueSequenceList u = uniqueSequenceList.get(0);
		Atom[] referenceAtoms = u.getReferenceChain();		
		AFPChain afp = alignPairBySequence(referenceAtoms, cAlphaAtoms);
		double identity = 0;
		double rmsd = 0;
		if (afp != null) {
			int[][][] al = afp.getOptAln();
			if (al != null) {
				identity = afp.getIdentity();
				rmsd = afp.getChainRmsd();
			}
		} else {
			System.out.println("AFPChain is null");
		}

		System.out.println("Smith-Waterman: seq. identity: " + (float)identity + " RMSD: " + (float)rmsd + " alignment length: " + afp.getOptLength());
//		System.out.println(afp.getAlnseq1());
//		System.out.println(afp.getAlnseq2());
		if (identity < 1.0) {
			// when identity is less than 1, there may be miss-alignments due to gaps (missing residues).
			// In that case use structural alignment.
			afp = alignPairByStructure(referenceAtoms, cAlphaAtoms);
			if (afp == null) {
				return false;
			}
			identity = afp.getIdentity();
			rmsd = afp.getChainRmsd();
			System.out.println("CE            : seq. identity: " + (float)identity + " RMSD: " + (float)rmsd + " alignment length: " + afp.getOptLength());
//			System.out.println(afp.getAlnseq1());
//			System.out.println(afp.getAlnseq2());
		}
		int[][][] alignment = afp.getOptAln();
		if (alignment != null) {
			UniqueSequenceList seqList = new UniqueSequenceList(cAlphaAtoms, chainId, modelNumber, sequence);
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
//			System.out.println("Sequence alignment: ");
//			System.out.println(align1);
//			System.out.println(align2);
			addUniqueSequenceList(seqList);
			return true;
		}
		return false;
	}

	private AFPChain alignPairBySequence(Atom[] ca1Seq, Atom[] ca2Seq) {
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
	
	private AFPChain alignPairByStructure(Atom[] ca1Seq, Atom[] ca2Seq) {
	    CeParameters params = new CeParameters();
		params.setMaxGapSize(-1);
		// should set this only when seq. id. is high
		params.setScoringStrategy(CeParameters.SEQUENCE_CONSERVATION);
		params.setSeqWeight(2.0);
		
        AFPChain afp = null;
		try {
			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			afp = algorithm.align(ca1Seq,ca2Seq,params);
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}            
		return afp;
	}
	
	
	private int alignByAtomSequence(Atom[] ca1Seq, Atom[] ca2Seq, List<Integer> align1, List<Integer> align2) {
		AFPChain afp = alignPairBySequence(ca1Seq, ca2Seq);
		int[][][] align = afp.getOptAln();
		if (align == null) {
			return 0;
		}
		int len =  afp.getOptLength();
//		System.out.println("identity: " + afp.getIdentity() + " len: " + len);

		List<Integer> delta = new ArrayList<Integer>();
		Set<Integer> unique = new HashSet<Integer>();

		for (int i = 0; i < len; i++) {
			Atom a1 = ca1Seq[align[0][0][i]];
			String residueName1 = a1.getGroup().getPDBName();
			Atom a2 = ca2Seq[align[0][1][i]];
			String residueName2 = a2.getGroup().getPDBName();
			if (residueName1.equals(residueName2)) {
			    int n1 = a1.getGroup().getResidueNumber().getSeqNum();
			    int n2 = a2.getGroup().getResidueNumber().getSeqNum();
			    delta.add(n2-n1);
			    unique.add(n2-n1);
			}
		}
		
		int offset = 0;
		int frequency = 0;
        for (Integer i: unique) {
        	int freq = Collections.frequency(delta, i);
        	if (freq > frequency) {
        		offset = i;
        		frequency = freq;
        	}
        }
        
        for (int i = 0; i < len; i++) {
        	Atom a1 = ca1Seq[align[0][0][i]];
			int n1 = a1.getGroup().getResidueNumber().getSeqNum();
			Atom a2 = ca2Seq[align[0][1][i]];
			int n2 = a2.getGroup().getResidueNumber().getSeqNum();
			if (n2 - offset == n1) {
				align1.add(align[0][0][i]);
				align2.add(align[0][1][i]);
			}
        }
//        System.out.println("PDB alignment: ");
//        System.out.println(align1);
//        System.out.println(align2);
        return align1.size();
	}
	
	private void createAlignedCAlphaAtoms() {
		List<Integer> indices = getReferenceResidueIndices();
		alignmentLength = indices.size();
//		System.out.println("In common: " + indices);
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
//				System.out.println("Unaligned: " + unalignedAtoms.length);
//				System.out.println("Aligned: " + alignedAtoms.length);
				for (int j = 0; j < alignedAtoms.length; j++) {
					alignedAtoms[j] = unalignedAtoms[alignmentIndices.get(j)];
				}
//				System.out.println("Calpha: " + alignedAtoms.length);
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
//		alignedCBetaAtoms = new ArrayList<Atom[]>();
//		for (int j = 0; j < cbetas.size(); j++) {
//			cb = cbetas.get(j);
//			Atom[] cbatoms = new Atom[inCommon.size()];
//			int index = 0;
//			for (int i = 0; i < cb.length; i++) {
//			   if (inCommon.contains(i)) {
////				   if (cb[i] == null) {
////					   System.out.println("Adding null: " + index);
////				   }
//				   cbatoms[index] = cb[i];
//				   index++;
//			   }
//			}
////			System.out.println("index: " + index);
//			alignedCBetaAtoms.add(cbatoms);
//			System.out.println("Cbeta: " + cbatoms.length);
//		}
	}
	
	private List<Integer> getReferenceResidueIndices() {
		List<Integer> indices = new ArrayList<Integer>(uniqueSequenceList.get(0).getAlignment1());
		for (UniqueSequenceList u: uniqueSequenceList) {
           indices.retainAll(u.getAlignment1());
		}
		return indices;
	}
}
