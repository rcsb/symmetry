package org.biojava3.structure.align.symm.quaternary;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;

public class UniqueSequenceList {
	private String sequenceString = "";
    private List<Integer> alignment1 = null;
    private List<Integer> alignment2 = null;
    private List<Atom[]> caAtoms = new ArrayList<Atom[]>();
    private List<String> chainIds = new ArrayList<String>();
    
    public UniqueSequenceList(Atom[] cAlphaAtoms, String chainId) {
    	this.caAtoms.add(cAlphaAtoms);
    	this.chainIds.add(chainId);
    	this.sequenceString =  getSequenceString(cAlphaAtoms);
    	this.alignment1 = new ArrayList<Integer>(cAlphaAtoms.length);
    	for (int i = 0; i < cAlphaAtoms.length; i++) {
    		this.alignment1.add(i);
    	}
    	this.alignment2 = alignment1;
    }
    
    /**
     * Return true is the sequence and residues numbers of the passed in array of
     * atoms matches those of this unique sequence list
     * 
     * @param caAlphaAtoms
     * @return
     */
    public boolean isMatch(Atom[] caAlphaAtoms) {
    	return sequenceString.equals(getSequenceString(caAlphaAtoms));
    }
    
    public boolean addChain(Atom[] cAlphaAtoms, String chainId) {
    	if (isMatch(cAlphaAtoms)) {
    		this.caAtoms.add(cAlphaAtoms);
        	this.chainIds.add(chainId); 
    		return true;
    	} else {
    		return false;
    	}
    }
    
    public int getChainCount() {
    	return caAtoms.size();
    }
    
    public Atom[] getChain(int index) {
    	return caAtoms.get(index);
    }
    
    public String getChainId(int index) {
    	return chainIds.get(index);
    }
    
    public Atom[] getReferenceChain() {
    	return caAtoms.get(0);
    }
    
     /**
	 * @return the sequenceString
	 */
	public String getSequenceString() {
		return sequenceString;
	}
	/**
	 * @param sequenceString the sequenceString to set
	 */
	public void setSequenceString(String sequenceString) {
		this.sequenceString = sequenceString;
	}
	/**
	 * @return the alignment1
	 */
	public List<Integer> getAlignment1() {
		return alignment1;
	}
	/**
	 * @param alignment1 the alignment1 to set
	 */
	public void setAlignment1(List<Integer> alignment1) {
		this.alignment1 = alignment1;
	}
	/**
	 * @return the alignment2
	 */
	public List<Integer> getAlignment2() {
		return alignment2;
	}
	/**
	 * @param alignment2 the alignment2 to set
	 */
	public void setAlignment2(List<Integer> alignment2) {
		this.alignment2 = alignment2;
	}
	
	public static String getSequenceString(Atom[] caAlphaAtoms) {
		StringBuilder builder = new StringBuilder();

		for (Atom a:  caAlphaAtoms) {
			Group g = a.getGroup();
			if (! g.getPDBName().equals("UNK")) {
				builder.append(g.getResidueNumber());
				builder.append(g.getPDBName());
			}
		}
		
		return builder.toString();
	}
     
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("#: ");
		builder.append(caAtoms.size());
		builder.append(" seq: ");
		builder.append(sequenceString);
		builder.append("\n");
		builder.append(alignment1);
		builder.append("\n");
		builder.append(alignment2);
		return builder.toString();
	}
	
}
