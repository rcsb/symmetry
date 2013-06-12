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
import org.biojava.bio.structure.align.util.AFPChainScorer;

public class SequenceAlignmentCluster {
	private QuatSymmetryParameters parameters = null;
	private List<UniqueSequenceList> uniqueSequenceList = new ArrayList<UniqueSequenceList>();
	private List<Atom[]> alignedCAlphaAtoms = null;

	private int alignmentLength = 0;
	private boolean pseudoStoichiometric = false;
	private double minSequenceIdentity = 1.0;
	private double maxSequenceIdentity = 0.0;
	
	private boolean modified = true;

	public SequenceAlignmentCluster (QuatSymmetryParameters parameters) {
		this.parameters = parameters;
	}
	
	public boolean isPseudoStoichiometric() {
		return pseudoStoichiometric;
	}
	
	public double getMinSequenceIdentity() {
		if (! isPseudoStoichiometric()) {
			return 1.0;
		}
		return minSequenceIdentity;
	}

	public double getMaxSequenceIdentity() {
		if (! isPseudoStoichiometric()) {
			return 1.0;
		}
		return maxSequenceIdentity;
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
//			System.out.println("in common: "  + inCommon);

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

	public int[][][] alignClustersByStructure(SequenceAlignmentCluster cluster) {
		// create aligned C alpha atom lists
		run();
		cluster.run();
		pseudoStoichiometric = false;
		
		int[][][] alignment = null;
		Atom[] referenceAtoms1 = this.getUniqueSequenceList().get(0).getCalphaAtoms();
		Atom[] referenceAtoms2 = cluster.getUniqueSequenceList().get(0).getCalphaAtoms();
		
		double alignmentLengthFraction = (double)Math.min(referenceAtoms1.length, referenceAtoms2.length) /
				Math.max(referenceAtoms1.length, referenceAtoms2.length);
		
		if (parameters.isVerbose()) {
			System.out.println("SequenceAlignmentCluster: alignmentLengthFraction: " + alignmentLengthFraction);
		}
		if (alignmentLengthFraction < parameters.getAlignmentFractionThreshold()) {
			return alignment;
		}
		
		AFPChain afp = alignPairByStructure(referenceAtoms1, referenceAtoms2);
		if (afp == null) {
			return alignment;
		}
		double identity = afp.getIdentity();
		double rmsd = afp.getChainRmsd();
		
		if (parameters.isVerbose()) {
			System.out.println("SequenceAlignmentCluster: CE: seq. identity: " + (float)identity + " RMSD: " + (float)rmsd + " alignment length: " + afp.getOptLength());
		}
		
		
	
		
//		System.out.println(afp.getAlnseq1());
//		System.out.println(afp.getAlnseq2());
	
    	alignment = afp.getOptAln();
    	if (alignment != null) {		
    		alignmentLengthFraction = (double)afp.getOptLength()/Math.max(referenceAtoms1.length, referenceAtoms2.length);
    		if (parameters.isVerbose()) {
    			System.out.println("SequenceAlignmentCluster: alignmentLengthFraction: " + alignmentLengthFraction);
    		}
//    		if (rmsd > parameters.getRmsdThreshold() || 
//    				alignmentLengthFraction < parameters.getAlignmentFractionThreshold()) {
//    			alignment = null;
//    			return alignment;
//    		}
    		
    		// alternative: tmSCore:
    		// double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
    		// if ( tmScore < 0.35) {
    		// return null ...
    		//}
    		
    		if ( !afp.isSignificantResult()  || 
    				alignmentLengthFraction < parameters.getAlignmentFractionThreshold() ) {
    			alignment = null;
    			return null;
    		}
    		
    	}
		
		List<Integer> align1 = new ArrayList<Integer>(alignment[0][0].length);
		for (Integer a1: alignment[0][0]) {
			align1.add(a1);
		}

		List<Integer> align2 = new ArrayList<Integer>(alignment[0][1].length);
		for (Integer a2: alignment[0][1]) {
			align2.add(a2);
		}

		// if pseudo stoichiometric, keep track of min and max idenity found
		if (identity < parameters.getSequencePseudoSymmetryThreshold()) {
			if (identity < minSequenceIdentity) {
				minSequenceIdentity = identity;
			}
			if (identity > maxSequenceIdentity) {
				maxSequenceIdentity = identity;
			}
			pseudoStoichiometric = true;
		}
		return alignment;
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
		//params.setMaxGapSize(-1);
		// should set this only when seq. id. is high
		//params.setScoringStrategy(CeParameters.SEQUENCE_CONSERVATION);
		//params.setSeqWeight(2.0);
		
        AFPChain afp = null;
		try {
			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			afp = algorithm.align(ca1Seq,ca2Seq,params);
			//System.out.println(afp.toFatcat(ca1Seq, ca2Seq));
		} catch (StructureException e) {
			// TODO Auto-generated catch block
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
