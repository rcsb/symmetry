package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;

public class ProteinChainExtractor  {
	private Structure structure = null;
	private int minSequenceLength = 24;
	private boolean unknownSequence = false;
	private boolean modified = true;
	
	private List<Atom[]> cAlphaTrace = new ArrayList<Atom[]>();	
	private List<String> chainIds = new ArrayList<String>();
	private List<Integer> modelNumbers = new ArrayList<Integer>();
	private List<String> sequences = new ArrayList<String>();

	public ProteinChainExtractor(Structure structure, int minSequenceLength) {
		this.structure = structure;
		this.minSequenceLength = minSequenceLength;
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
	
	public boolean containsUnknownSequence() {
		return unknownSequence;
	}
	
    private void run() {
    	if (modified) {
    		extractProteinChains();
    		modified = false;
    	}
    }

	private void extractProteinChains() {
		int models = 1;
		if (structure.isBiologicalAssembly()) {
			models = structure.nrModels();
		}
		
		for (int i = 0; i < models; i++) {
			for (Chain c : structure.getChains(i)) {
				Atom[] ca = StructureTools.getAtomCAArray(c);
				if (containsUnknownResidues(ca)) {
					ca = removeUnknownResidues(ca);
					unknownSequence = true;
				}
				if (ca.length >= minSequenceLength) {
				   System.out.println("Chain " + c.getChainID() + ": " + c.getSeqResSequence());
//				   System.out.println("Number CA atoms: " + ca.length);
				   cAlphaTrace.add(ca);
				   chainIds.add(c.getChainID());
				   modelNumbers.add(i);
				   sequences.add(c.getSeqResSequence());
				}
			}
		}
	}
	
	private boolean containsUnknownResidues(Atom[] atoms) {
		for (Atom atom: atoms) {
			if (atom.getGroup().getPDBName().equalsIgnoreCase("UNK")) {
				return true;
			}
		}
		return false;
	}
	
	private Atom[] removeUnknownResidues(Atom[] atoms) {
		List<Atom> atomList = new ArrayList<Atom>(atoms.length);
		for (Atom atom: atoms) {
			if (! atom.getGroup().getPDBName().equalsIgnoreCase("UNK")) {
				atomList.add(atom);
			}
		}
		return atomList.toArray(new Atom[atomList.size()]);
	}
}
