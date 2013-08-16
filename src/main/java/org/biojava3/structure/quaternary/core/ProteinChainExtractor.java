package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;

public class ProteinChainExtractor  {
	private Structure structure = null;
	private QuatSymmetryParameters parameters = null;
	private boolean modified = true;
	private int adjustedMinimumSequenceLength = 0;
	
	private List<Atom[]> cAlphaTrace = new ArrayList<Atom[]>();	
	private List<String> chainIds = new ArrayList<String>();
	private List<Integer> modelNumbers = new ArrayList<Integer>();
	private List<String> sequences = new ArrayList<String>();

	public ProteinChainExtractor(Structure structure, QuatSymmetryParameters parameters) {
		this.structure = structure;
		this.parameters = parameters;
		modified = true;
	}
	
	public List<Atom[]> getCalphaTraces() {
		run();
		return cAlphaTrace;
	}
	
	public List<String> getChainIds() {
        run();
		return chainIds;
	}
	
	public List<Integer> getModelNumbers() {
        run();
		return modelNumbers;
	}
	
	public List<String> getSequences() {
		run();
		return sequences;
	}
	
	public int getAdjustedMinimumSequenceLength() {
		run();
		return adjustedMinimumSequenceLength;
	}
	
    private void run() {
    	if (modified) {
    		calcAdjustedMinimumSequenceLength();
    		extractProteinChains();
    		modified = false;
    	}
    }

	private void extractProteinChains() {
		int models = 1;
		if (structure.isBiologicalAssembly()) {
			models = structure.nrModels();
		}
		
		if (parameters.isVerbose()) {
			System.out.println("Protein chains used in calculation:");
			System.out.println("Adjusted minimum sequence length: " + adjustedMinimumSequenceLength);
		}
		
		for (int i = 0; i < models; i++) {
			for (Chain c : structure.getChains(i)) {
				Atom[] ca = StructureTools.getAtomCAArray(c);
				ca = retainStandardAminoAcidResidues(ca);

				if (ca.length >= adjustedMinimumSequenceLength) {
					if (parameters.isVerbose()) {
				        System.out.println("Chain " + c.getChainID() + " Calpha atoms: " + ca.length + " seqres: " + c.getSeqResSequence());
					}
				  
				   cAlphaTrace.add(ca);
				   chainIds.add(c.getChainID());
				   modelNumbers.add(i);
				   sequences.add(replaceQuestionMarks(c.getSeqResSequence()));
				}
			}
		}
	}
	
	/**
	 * Returns an adapted minimum sequence length. This method ensure that
	 * structure that only have short chains are not excluded by the
	 * minimumSequenceLength cutoff value;
	 * @return
	 */
	private void calcAdjustedMinimumSequenceLength() {
		int models = 1;
		if (structure.isBiologicalAssembly()) {
			models = structure.nrModels();
		}
		
		int maxLength = Integer.MIN_VALUE;
		int minLength = Integer.MAX_VALUE;
		
		List<Integer> lengths = new ArrayList<Integer>();
		for (int i = 0; i < models; i++) {
			for (Chain c : structure.getChains(i)) {
				Atom[] ca = StructureTools.getAtomCAArray(c);
				ca = retainStandardAminoAcidResidues(ca);
				if (ca.length >= parameters.getAbsoluteMinimumSequenceLength()) {
					maxLength = Math.max(ca.length, maxLength);
					minLength = Math.min(ca.length, minLength);
					lengths.add(ca.length);
				}
			}
		}
		
		
		
		adjustedMinimumSequenceLength = parameters.getMinimumSequenceLength();
		
		if (lengths.size() < 2) {
			return;
		}
		
		double median = 0;
		Collections.sort(lengths);
		if (lengths.size() %2 == 1) {
			int middle = (lengths.size()-1) / 2;
			median = lengths.get(middle);
		} else {
	        int middle2 = lengths.size()/2;
	        int middle1 = middle2-1;
			median = 0.5 * (lengths.get(middle1) + lengths.get(middle2));
		}
	
		if (minLength >= median * parameters.getMinimumSequenceLengthFraction()) {
			adjustedMinimumSequenceLength = Math.min(minLength, parameters.getMinimumSequenceLength());
		}
	}
	
	// In some cases "?" are in the sequence string. Example 2WS1. This is caused
	// because the chemical component file YNM doesn't contain a one-letter code.
	private String replaceQuestionMarks(String sequence) {
		return sequence.replaceAll("\\?", "X");
	}
	
	private Atom[] retainStandardAminoAcidResidues(Atom[] atoms) {
		List<Atom> atomList = new ArrayList<Atom>(atoms.length);
		for (Atom atom: atoms) {
			Group group = atom.getGroup();
			if (group.getPDBName().equalsIgnoreCase("UNK")) {
				continue;
			}
			if (! group.getType().equals("amino") && ! isSpecialAminoAcid(group)) {
				System.out.println("Excluding: " + group.getPDBName());
				continue;
			}
			atomList.add(atom);
		}
		return atomList.toArray(new Atom[atomList.size()]);
	}
		
	private boolean isSpecialAminoAcid(Group group) {
		String groupName = group.getPDBName();
		return groupName.equals("MSE") || groupName.equals("HYP") || groupName.equals("HZP") ||
				groupName.equals("MP8") || groupName.equals("FP9") || groupName.equals("DPR") || groupName.equals("NLE");
	}
}
