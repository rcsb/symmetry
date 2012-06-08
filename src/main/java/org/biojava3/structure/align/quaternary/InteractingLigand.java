package org.biojava3.structure.align.quaternary;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;

public class InteractingLigand {
     private Group ligand = null;
     private List<Chain> chains = new ArrayList<Chain>();
     private List<String> compositionIds = new ArrayList<String>();
     private List<Integer> contacts = new ArrayList<Integer>();
     
     public InteractingLigand(Group ligand) {
    	 this.ligand = ligand;
     }
     
     public void addInteraction(Chain chain, String compositionId, int interactions) {
    	 System.out.println("Adding interactions: " + chain.getChainID() + " compositionId: " + compositionId + " contacts: " + interactions);
    	 chains.add(chain);
    	 compositionIds.add(compositionId);
    	 contacts.add(interactions);
     }
     
     public String getLigandCode() {
    	 return ligand.getPDBName();
     }
     
     public int getInteractingSubunitCount() {   
    	return chains.size();
     }
     
     public int getTotalInteractions() {
    	 int total = 0;
    	 for (Integer c: contacts) {
    		 total += c;
    	 }
    	 return total;
     }
     
     public String chainIdSignature() {
    	 StringBuilder builder = new StringBuilder();
    	 builder.append(ligand.getPDBName());
    	 builder.append("[");
    	 for (int i = 0; i < chains.size(); i++) {
    		 builder.append(chains.get(i).getChainID());

    	 }
    	 builder.append("]");
    	 return builder.toString();
     }
     
     public String toString() {
    	 StringBuilder builder = new StringBuilder();
    	 builder.append(ligand.getPDBName());
    	 builder.append("[");
    	 for (int i = 0; i < chains.size(); i++) {
    		 builder.append(compositionIds.get(i));

    	 }
    	 builder.append("]");
    	 return builder.toString().trim();
     }
}
