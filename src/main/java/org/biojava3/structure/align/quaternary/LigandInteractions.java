package org.biojava3.structure.align.quaternary;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Compound;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava3.structure.align.symm.quaternary.analysis.ProteinComplexSignature;

public class LigandInteractions {
	private Structure structure = null;
	private ProteinComplexSignature signature = null;
	private List<Chain> chains = new ArrayList<Chain>(0);
	private List<InteractingLigand> ligands = new ArrayList<InteractingLigand>();
	private boolean modified = true;
	private int interactionType = 0;

	public LigandInteractions(Structure structure, ProteinComplexSignature signature) {
		this.structure = structure;
		this.signature = signature;
		modified = true;
	}
	
	public void setInteractingChains(List<Chain> chains) {
		this.chains = chains;
		modified = true;
		System.out.println("Number of chains: " + chains.size());
	}
	
	public int getInteractionType() {
		run();
		return interactionType;
	}
	
//	public List<InteractingLigand> getInteractingLigands() {
//		run();
//		return ligands;
//	}
	
	private List<Group> getLigands() {
		List<Group> ligs = new ArrayList<Group>();
		// TODO ligands could be in their own chains, they would be ignored here!!!
		for (Chain c: structure.getChains()) {
			for (Group g: c.getAtomLigands()) {
				String type = g.getType();
				double mass = getGroupMass(g);
				if (type.equals("hetatm") && mass > 100) {
					ligs.add(g);
					System.out.println("Ligand: " + g + " mass: " + mass);
				}
			}
		}	
		return ligs;
	}
	
	private void run() {	
		if (modified) {
			interactionType = 0;
			boolean intraSubunit = false;
			boolean interSubunit = false;
			for (Group l: getLigands()) {
				// calculate number of contact that each ligand makes with protein chains
				InteractingLigand lig = new InteractingLigand(l);
				//			System.out.println("Ligand: " + l);
				for (Chain c: chains) {
					int contacts = getContacts(c,l);
					if (contacts > 5) {
						String compositionId = signature.getCompositionId(c.getChainID());
						lig.addInteraction(c, compositionId, contacts);
					}
				}

				int count = lig.getInteractingSubunitCount();
				if (count == 1) {
					intraSubunit = true;
				} else if (count > 1) {
					interSubunit = true;
				}
				
				if (count > 0) {
					ligands.add(lig);
				}
			}
			if (intraSubunit) {
				interactionType += 1;
			}
			if (interSubunit) {
				interactionType += 2;
			}
			modified = false;
		}
	}
	   
	public String toString() {
		run();
		StringBuilder builder = new StringBuilder();
		Map<String,Integer> map = createSignatureMap();
		for (Entry<String, Integer> entry: map.entrySet()) {
			int value = entry.getValue();
			if (value > 1) {
				builder.append("(");
			}
			builder.append(entry.getKey());
			if (value > 1) {
				builder.append(")");
				builder.append(value);
			}
			builder.append(" ");
		}
		return builder.toString();
	}
	
	private Map<String,Integer> createSignatureMap() {
		Map<String,Integer> map = new TreeMap<String,Integer>();
		for (InteractingLigand lig: ligands) {
			Integer count = map.get(lig.toString());
			if (count == null) {
				count = 1;
			} else {
				count++;
			}
			map.put(lig.toString(),count);
		}
		return map;
	}

	private int getContacts(Chain chain, Group ligand) {
		int contacts = 0;
		for (Atom l: ligand.getAtoms()) {
			boolean inContact = false;
			for (Group g: chain.getAtomGroups()) {
				for (Atom a: g.getAtoms()) {
					try {
						if (Calc.getDistanceFast(a, l) < 25) {
							inContact = true;
							break;
						}
					} catch (StructureException e) {
					}
				}
				if (inContact) {
					break;
				}
			}
			if (inContact) {
				contacts++;
			}
		}
		return contacts;
	}
	
	private double getGroupMass(Group group) {
		double mass = 0.0f;
		for (Atom a: group.getAtoms()) {
			mass += a.getElement().getAtomicMass();
		}
		return mass;
	}
}
