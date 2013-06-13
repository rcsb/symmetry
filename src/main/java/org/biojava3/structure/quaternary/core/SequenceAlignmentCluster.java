package org.biojava3.structure.quaternary.core;

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
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;

public class SequenceAlignmentCluster implements Cloneable {
	private QuatSymmetryParameters parameters = null;
	private List<UniqueSequenceList> uniqueSequenceList = new ArrayList<UniqueSequenceList>();
	private List<Atom[]> alignedCAlphaAtoms = null;

	private int alignmentLength = 0;
	private double minSequenceIdentity = 1.0;
	private double maxSequenceIdentity = 0.0;
	
	private boolean modified = true;

	public SequenceAlignmentCluster (QuatSymmetryParameters parameters) {
		this.parameters = parameters;
	}
	
	public boolean isPseudoStoichiometric() {
		return minSequenceIdentity < parameters.getSequencePseudoSymmetryThreshold();
	}
	
	public double getMinSequenceIdentity() {
		if (! isPseudoStoichiometric()) {
			return 1.0;
		}
		return minSequenceIdentity;
	}
	
	public void setMinSequenceIdentity(double minSequenceIdentity) {
		this.minSequenceIdentity = minSequenceIdentity;
	}

	public double getMaxSequenceIdentity() {
		if (! isPseudoStoichiometric()) {
			return 1.0;
		}
		return maxSequenceIdentity;
	}

	public void setMaxSequenceIdentity(double maxSequenceIdentity) {
		this.maxSequenceIdentity = maxSequenceIdentity;
	}
	
	public void addUniqueSequenceList(UniqueSequenceList sequenceList) {
		uniqueSequenceList.add(sequenceList);
		modified = true;
	}

	public int getSequenceCount() {
		return uniqueSequenceList.size();
	}
	
	public int getSequenceAlignmentLength() {
		run();
		return alignmentLength;
	}
	
	public List<UniqueSequenceList> getUniqueSequenceList() {
		return uniqueSequenceList;
	}

	public List<String> getChainIds() {
		List<String> ids = new ArrayList<String>();
		for (UniqueSequenceList list: uniqueSequenceList) {
			ids.add(list.getChainId());
		}
		return ids;
	}
	
	public List<Integer> getModelNumbers() {
		List<Integer> numbers = new ArrayList<Integer>();
		for (UniqueSequenceList list: uniqueSequenceList) {
			numbers.add(list.getModelNumber());
		}
		return numbers;
	}
	
	public List<Integer> getStructureIds() {
		List<Integer> numbers = new ArrayList<Integer>();
		for (UniqueSequenceList list: uniqueSequenceList) {
			numbers.add(list.getStructureId());
		}
		return numbers;
	}
	
	public List<Atom[]> getAlignedCalphaAtoms() {
		run();
		return alignedCAlphaAtoms;
	}
	
	public boolean identityMatch(Atom[] cAlphaAtoms, String chainId, int modelNumber, int structureId, String sequence) {
		UniqueSequenceList u = uniqueSequenceList.get(0);

		// check for 100% identity match of reference sequence
		String refSequence = u.getSeqResSequence();
		boolean seqMatch = refSequence.equals(sequence);

		// if reference (SEQRES) sequences match 100%, 
		// find alignment of atom sequences by Smith-Waterman alignment
		if (seqMatch) {
			List<Integer> alig1 = new ArrayList<Integer>();
			List<Integer> alig2 = new ArrayList<Integer>();
			Atom[] referenceAtoms = u.getCalphaAtoms();
			int inCommon = alignIdenticalSequence(referenceAtoms, cAlphaAtoms, alig1, alig2);

			if (inCommon > 0) {
				UniqueSequenceList seqList = new UniqueSequenceList(cAlphaAtoms, chainId, modelNumber, structureId, sequence);
				seqList.setAlignment1(alig1);
				seqList.setAlignment2(alig2);
				//			System.out.println(alig1);
				//			System.out.println(alig2);
				addUniqueSequenceList(seqList);
				return true;
			}
		}

		return false;
	}
	
	public PairwiseAlignment getPairwiseAlignment(SequenceAlignmentCluster other) {
		PairwiseAlignment alignment = new PairwiseAlignment(this, other);
		
		Atom[] referenceAtoms1 = this.getUniqueSequenceList().get(0).getCalphaAtoms();
		Atom[] referenceAtoms2 = other.getUniqueSequenceList().get(0).getCalphaAtoms();
		
		double alignmentLengthFraction = (double)Math.min(referenceAtoms1.length, referenceAtoms2.length) /
				Math.max(referenceAtoms1.length, referenceAtoms2.length);
	
		if (alignmentLengthFraction < parameters.getAlignmentFractionThreshold()) {
			return null;
		}
		
		AFPChain afp = alignPairByStructure(referenceAtoms1, referenceAtoms2);
		if (afp == null) {
			return null;
		}
		
		if (! afp.isSignificantResult()) {
			return null;

    		// alternative: tmSCore:
    		// double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
    		// if ( tmScore < 0.35) {
    		// return null ...
		}
		
		int[][][] align = afp.getOptAln();
		if (align == null) {
			return null;
		}
			
    	alignmentLengthFraction = (double)afp.getOptLength()/Math.max(referenceAtoms1.length, referenceAtoms2.length);
    	if (parameters.isVerbose()) {
    			System.out.println("SequenceAlignmentCluster: alignmentLengthFraction: " + alignmentLengthFraction);
    	}
    	alignment.setAlignmentLengthFraction(alignmentLengthFraction);
    	alignment.setRmsd(afp.getChainRmsd());
    	alignment.setSequenceIdentity(afp.getIdentity());
    	alignment.setAlignment(afp.getOptAln());
    	
		return alignment;
	}
	
	public Object clone() {
	    SequenceAlignmentCluster copy = null;
		try {
			copy = (SequenceAlignmentCluster) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		// deep copy sequences
		copy.uniqueSequenceList = new ArrayList<UniqueSequenceList>();
		for (UniqueSequenceList seq: this.getUniqueSequenceList()) {
			copy.addUniqueSequenceList((UniqueSequenceList) seq.clone());
		}
		return copy;
	}
	
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
			alignedCAlphaAtoms = null;
			createAlignedCAlphaAtoms();
			modified = false;
		}
	}

	private AFPChain alignPairBySequence(Atom[] ca1Seq, Atom[] ca2Seq) {
		SmithWaterman3Daligner aligner = new SmithWaterman3Daligner();
		AFPChain afp = null;
		try {
			afp = aligner.align(ca1Seq, ca2Seq);
		} catch (StructureException e) {
			e.printStackTrace();
		} 
		return afp;
	}
	
	private AFPChain alignPairByStructure(Atom[] ca1Seq, Atom[] ca2Seq) {
       CeParameters params = new CeParameters();

        AFPChain afp = null;
		try {
			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			afp = algorithm.align(ca1Seq,ca2Seq,params);
			if (parameters.isVerbose()) {
				System.out.println(afp.toFatcat(ca1Seq, ca2Seq));
			}
		} catch (StructureException e) {
			e.printStackTrace();
		}            
		return afp;
	}
	
	
	private int alignIdenticalSequence(Atom[] ca1Seq, Atom[] ca2Seq, List<Integer> align1, List<Integer> align2) {
		AFPChain afp = alignPairBySequence(ca1Seq, ca2Seq);
		int[][][] align = afp.getOptAln();
		if (align == null) {
			return 0;
		}
		int len =  afp.getOptLength();
		if (parameters.isVerbose()) {
			System.out.println("SequenceAlignmentCluster: Smith-Waterman alignment: seq. identity:  " + afp.getIdentity());
		}
//		double identity = afp.getIdentity();
//		setMinSequenceIdentity(Math.min(getMinSequenceIdentity(),  identity));
//		setMaxSequenceIdentity(Math.max(getMaxSequenceIdentity(),  identity));

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
			Atom[] unalignedAtoms = u.getCalphaAtoms();
			Atom[] alignedAtoms = new Atom[alignmentIndices.size()];
			for (int j = 0; j < alignedAtoms.length; j++) {
				alignedAtoms[j] = unalignedAtoms[alignmentIndices.get(j)];
			}
			alignedCAlphaAtoms.add(alignedAtoms);
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
