package org.biojava.nbio.structure.align.symm.subunit;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;

public class MultipleAFP {

	List<AFPChain> allAlignments;
	List<List<Integer>> groups;
	Atom[] ca1;
	
	//List to store the residues aligned in the blocks. Dimensions are: [order][block_number][length of the block]
	List<List<List<Integer>>> blocks;
	//List to store the residues that are part of the free pool. Dimensions are: [order][pool_number][length of the pool]
	List<List<Integer>> free_pool;
	
	//This class pretends to use the CEMC approach for multiple structural alignment using the pairwise alignments obtained
	//from the align method.
	//TODO
	
	
	/**
	 * Method that saves a multiple alignment of the subunits in a file in FASTA format.
	 * 
	 * INPUT: the groups of residues sorted increasingly and the Aom[] array of the protein.
	 */
	public static void saveMultipleAln(List<List<Integer>> groups, Atom[] ca1) throws StructureException, IOException{
		
		int order = groups.get(0).size();
		String sFileName = "/home/scratch/align/unknown.fasta";
		
	    FileWriter writer = new FileWriter(sFileName);
		
		String[] subunits = new String[order];  //the alignment string of every subunit, with its gaps
		Arrays.fill(subunits, "");
		int[] position = new int[order];  //the positions in every subunit
		int[] next_position = new int[order];
		Arrays.fill(position, 0);
		for (int j=0; j<order; j++){
			next_position[j] = groups.get(position[j]).get(j);
		}
		
	    //Loop for every residue to see if it is included or not in the alignments and if there is a gap
		while (true){
			boolean stop = false;
			int gaps = 0;
			char[] provisional = new char[order];
			
			//If the position is higher than the subunit insert a gap
			for (int j=0; j<order; j++){
				if (position[j]>groups.size()-1){
					provisional[j] = '-';
					gaps++;
				}
				else {
					//If the next position is lower than the residue aligned there is a gap, so increment the gap
					int res = groups.get(position[j]).get(j);
					if (next_position[j]<res){
						provisional[j] = StructureTools.get1LetterCode(ca1[next_position[j]].getGroup().getPDBName());
						gaps++;
					}
					//If they are the same do not increment gap and consider a gap in case other subunits have residues in between
					else {
						provisional[j] = '-';
					}
				}
			}
			//If all sequences have gaps means that there are unaligned residues, so include them in the alignment
			if (gaps==order){
				for (int j=0; j<order; j++){
					subunits[j] += StructureTools.get1LetterCode(ca1[next_position[j]].getGroup().getPDBName());
					next_position[j]++;
				}
			}
			//If there are not gaps add the aligned residues
			else if (gaps==0){
				for (int j=0; j<order; j++){
					subunits[j] += StructureTools.get1LetterCode(ca1[next_position[j]].getGroup().getPDBName());
					position[j]++;
					next_position[j]++;
				}
			}
			//If only some subunits have gaps consider this information and add gaps to the subunits with missing residues
			else{
				for (int j=0; j<order; j++){
					if (provisional[j] == '-'){
						subunits[j] += '-';
					}
					else{
						subunits[j] += StructureTools.get1LetterCode(ca1[next_position[j]].getGroup().getPDBName());
						next_position[j] ++;
					}
				}
			}
			//Stop if all of the subunits have been analyzed until the end (all residues in the group)
			stop = true;
			for (int q=0; q<order; q++){
				if (position[q] < groups.size())
					stop = false;
			}
			//Stop if any subunit has reached the end of the molecule
			for (int q=0; q<order; q++){
				if (next_position[q] > ca1.length-1)
					stop = true;
			}
			if (stop) break;
	    }
		//Store in a fasta file the alignment for further analysis
	    for (int k=0; k<order; k++){
	    	if (k==0){
	    		writer.append(">No_rotation\n");
	    		writer.append(subunits[k]+"\n");
	    		System.out.println(subunits[k]);
	    	}
	    	else{
	    		writer.append(">"+((360/order)*k)+"rotation\n");
	    		writer.append(subunits[k]+"\n");
	    		System.out.println(subunits[k]);
	    	}
	    }
	    writer.flush();
	    writer.close();
	}
}
