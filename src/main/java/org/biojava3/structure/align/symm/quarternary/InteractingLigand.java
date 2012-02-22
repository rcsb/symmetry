package org.biojava3.structure.align.symm.quarternary;

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
     
     public int getInteractingSubunitCount(float interactingFraction) {   
    	 int total = getTotalInteractions();
    	 System.out.println("total: " + total);
    	 if (total == 0) {
    		 return 0;
    	 }
    	 
    	 int subunitCount = 0;
    	 for (Integer c: contacts) {
    		 if ((float) c/(float)total > interactingFraction) {
    			 subunitCount++;
    		 }
    	 }
    	 System.out.println("Subunit count: " + subunitCount);
    	 return subunitCount;
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
    	 return builder.toString();
     }
}
