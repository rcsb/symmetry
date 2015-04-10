package org.biojava.nbio.structure.align.symm.subunit;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * 	This class pretends to use the CEMC approach for multiple structural alignment of the subunits
 *  using the refined pairwise alignment obtained from the align method. Optimization step.
 *  Also includes a method to obtain a multiple sequence alignment as text format.
 *  
 *  @author lafita
 *
 */
public class MultipleAFP {

	private static final boolean debug = true;
	
	AFPChain afpChain;
	Atom[] ca;
	int order;
	int subunitLen;
	
	//Variables that store the history of the optimization, in order to be able to plot the evolution of the system.
	List<Integer> subunitLenHistory;
	List<Double> rmsdHistory;
	List<Double> scoreHistory;
	
	//List to store the residues aligned, in the block. Dimensions are: [order][subunitLen]
	List<ArrayList<Integer>> block;
	//List to store the residues not aligned, that are part of the free pool. Dimensions are: [order][residues in the pool]
	List<ArrayList<Integer>> freePool;
	
	//They store the average RMSD and score of all possible rotation superpositions (alignments) of the protein
	double rmsd;
	double score;
	
	//Multiple alignment String sequences
	String[] alnSequences;
	String alnSymbols;
	
	public MultipleAFP(AFPChain afp, Atom[] ca1){
		
		//Initialize member variables
		afpChain = afp;
		ca = ca1;
		order = afpChain.getBlockNum();
		subunitLen = afpChain.getOptLen()[0];
		
		//Initialize alignment variables
		block = new ArrayList<ArrayList<Integer>>();
		freePool = new ArrayList<ArrayList<Integer>>();
		
		//Store the residues that have been added either to the block or to the freePool
		List<Integer> alreadySeen = new ArrayList<Integer>();
		
		//Generate the initial state of the system from the aligned blocks of the AFPChain
		for (int i=0; i<order; i++){
			ArrayList<Integer> residues = new ArrayList<Integer>();
			for (int j=0; j<subunitLen; j++){
				Integer residue = afpChain.getOptAln()[i][0][j];
				residues.add(residue);
				alreadySeen.add(residue);
			}
			block.add(residues);
			freePool.add(new ArrayList<Integer>());
		}
		
		//Add any residue not aligned to the free pool, in the corresponding subunit region (provisional initial state)
		for (int i=0; i<ca.length; i++){
			if (alreadySeen.contains(i)) continue;
			for (int j=0; j<order; j++){
				if (i<block.get(j).get(subunitLen-1) && !alreadySeen.contains(i)){
					freePool.get(j).add(i);
					alreadySeen.add(i);
				}
			}
			if (!alreadySeen.contains(i)){
				freePool.get(order-1).add(i);
			}
		}
		
		try {
			updateAvgRMSDandScore();
			//The maxIter is set in function of the protein length
			optimizeMC(100*ca.length);
			saveHistory("/scratch/mcopt/"+afpChain.getName1()+"_MCref.csv");
			if (debug) System.out.println("Saved history.");
			//saveSeqAln("/scratch/align/"+afpChain.getName1()+"_MC.fasta");
		} catch (StructureException e) {
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 *  Optimization method based in a Monte-Carlo approach. Starting from the refined afpChain uses 5 types of moves:
	 *  
	 *  	1- Shift Row: if there are enough freePool residues available.
	 *  	2- Change Residue: move residues from one subunit to another.
	 *  	3- Expand Block: if there are enough freePool residues available.
	 *  	4- Shrink Block: move a block column to the freePool.
	 *  	5- Split and Shrink Block: split a block in the middle and shrink one column.
	 */
	@SuppressWarnings("unchecked")
	private void optimizeMC(int maxIter) throws StructureException{
		
		//Initialize the history variables
		subunitLenHistory = new ArrayList<Integer>();
		rmsdHistory = new ArrayList<Double>();
		scoreHistory = new ArrayList<Double>();
		
		//Initialize a random generator number
		Random rnd = new Random();
		int conv = 0;  //Number of steps without an alignment improvement
		double C = 0.1;  //Constant for the probability of acceptance for bad moves, will be decreased
		
		int i = 1;
		
		while (i<maxIter && conv<(1000)){
			
			//Save the state of the system in case the modifications are not favorable
			List<ArrayList<Integer>> lastBlock = new ArrayList<ArrayList<Integer>>();
			List<ArrayList<Integer>> lastFreePool = new ArrayList<ArrayList<Integer>>();
			for (int k=0; k<order; k++){
				lastBlock.add((ArrayList<Integer>) block.get(k).clone());
				lastFreePool.add((ArrayList<Integer>) freePool.get(k).clone());
			}
			double lastRMSD = rmsd;
			double lastScore = score;
			
			moveResidue();  //At the beginning of each iteration move randomly a residue from the freePool
			
			//Randomly select one of the steps to modify the alignment
			int move = rnd.nextInt(4);
			switch (move){
			case 0: shiftRow();
					if (debug) System.out.println("did shift");
					break;
			case 1: expandBlock();
					if (debug) System.out.println("did expand");
					break;
			case 2: shrinkBlock();
					if (debug) System.out.println("did shrink");
					break;
			case 3: splitBlock();
					if (debug) System.out.println("did split");
					break;
			}
			
			//Get the properties of the new alignment
			updateAvgRMSDandScore();
			
			double AS = score-lastScore;  //Change in the alignment score
			double prob=0.0;
			
			if (AS<0){
				
				//Probability of accepting the new alignment
				prob = Math.min(Math.max((C+AS)/Math.pow(i,0.75),0.0),1.0);
				double p = rnd.nextDouble();
				//Reject the move
				if (p>prob){
					block = lastBlock;
					freePool = lastFreePool;
					rmsd = lastRMSD;
					score = lastScore;
					subunitLen = block.get(0).size();
					conv ++; //Increment the number of steps without a change in score
					
				} else conv = 0;
				
			} else if (AS!=0.0) conv=0;
			
			if (debug) 	System.out.println(i+": --prob: "+prob+", --score: "+AS+", --conv: "+conv);
			
			subunitLenHistory.add(subunitLen);
			rmsdHistory.add(rmsd);
			scoreHistory.add(score);
			
			i++;
		}
		
		int[][][] newAlgn = new int[order][2][subunitLen];
		for (int su=0; su<order; su++){
			int[] chain1 = new int[subunitLen];
			int[] chain2 = new int[subunitLen];
			for (int k=0; k<subunitLen; k++){
				chain1[k] = block.get(su).get(k);
				chain2[k] = block.get((su+1)%order).get(k);
			}
			newAlgn[su][0] = chain1;
			newAlgn[su][1] = chain2;
		}
		
		afpChain = AlignmentTools.replaceOptAln(newAlgn, afpChain, ca, ca);
		updateSeqAln();
	}
	
	/**
	 *  Move all the block residues of one subunit one position to the left or right and move the corresponding
	 *  boundary residues from the freePool to the block, and viceversa.
	 */
	private void shiftRow(){
		
		//Initialize a random generator number
		Random rnd = new Random();
		boolean shifted = false;
		int iter = 0;
		
		while(!shifted && iter<order){
			
			iter++;
			int su = rnd.nextInt(order); //Select randomly the subunit that is going to be shifted
			int rl = rnd.nextInt(2);  //Select between moving right (0) or left (1)
			int res = rnd.nextInt(subunitLen); //Residue as a pivot to make the shift
			
			if (freePool.get(su).size()==0) continue;  //If the freePool is empty the subunit cannot be shifted
			
			switch(rl){
			case 0: //Move to the right
				//Check that there is at least one residue in the freePool to the left (smaller) than the pivot
				if (freePool.get(su).get(0)<block.get(su).get(res)){
					
					int leftGap = res-1;  //Find the nearest gap to the left of the res
					for (int i = res; i>=0; i--){
						if(leftGap < 0){
							break;
						} else if (block.get(su).get(i) > block.get(su).get(leftGap)+1){
							break;
						}
						leftGap--;
					}
					
					int rightGap = res+1;  //Find the nearest gap to the right of the res
					for (int i = res; i<=subunitLen; i++){
						if(rightGap == subunitLen){
							break;
						} else if (block.get(su).get(i)+1 < block.get(su).get(rightGap)){
							break;
						}
						rightGap++;
					}
					
					//Move the residue at the left of the block from the freePool to the block
					Integer residue = block.get(su).get(leftGap+1)-1;
					block.get(su).add(leftGap+1,residue);
					freePool.get(su).remove(residue);
					
					//Move the residue at the right of the block to the freePool
					if (rightGap == subunitLen){
						freePool.get(su).add(block.get(su).get(rightGap));
						block.get(su).remove(rightGap);
					}
					else{
						freePool.get(su).add(block.get(su).get(rightGap+1));
						block.get(su).remove(rightGap+1);
					}
					Collections.sort(freePool.get(su));
					
					shifted = true;
				}
				break;
				
			case 1: //Move to the left
				//Check that there is at least one residue in the freePool to the right (bigger) than the pivot
				if (freePool.get(su).get(freePool.get(su).size()-1)>block.get(su).get(res)){
					
					int leftGap = res-1;  //Find the nearest gap to the left of the res
					for (int i = res; i>=0; i--){
						if(leftGap <= 0){
							leftGap = 0;
							break;
						} else if (block.get(su).get(i) > block.get(su).get(leftGap)+1){
							leftGap++;
							break;
						}
						leftGap--;
					}
					
					int rightGap = res+1;  //Find the nearest gap to the right of the res
					for (int i = res; i<=subunitLen; i++){
						if(rightGap >= subunitLen){
							rightGap = subunitLen;
							break;
						} else if (block.get(su).get(i)+1 < block.get(su).get(rightGap)){
							break;
						}
						rightGap++;
					}
					
					//Move the residue at the right of the block from the freePool to the block
					Integer residue = block.get(su).get(rightGap-1)+1;
					if (rightGap == subunitLen){
						block.get(su).add(residue);
						freePool.get(su).remove(residue);
					}
					else {
						block.get(su).add(rightGap,residue);
						freePool.get(su).remove(residue);
					}
					
					//Move the residue at the left of the block to the freePool
					freePool.get(su).add(block.get(su).get(leftGap));
					Collections.sort(freePool.get(su));
					block.get(su).remove(leftGap);
					
					shifted = true;
				}
				break;
			}
		}
	}
	
	/**
	 * Move a residue of the freePool from the start of one subunit to the end of the previous, or
	 *                                from the end of one subunit to the start of the next.
	 * Exclude the first and last residues, because makes no sense to move them.
	 * This move does not cause a score change, so it is made each iteration with combination to other moves.
	 */
	private void moveResidue(){
		
		//Initialize a random generator number
		Random rnd = new Random();
		int su = rnd.nextInt(order);
		int rl = rnd.nextInt(2);  //Randomly choose to move one residue from the right-end (0) or left-start (1)
		
		//Exclude first and last residues from being moved.
		if (su==0) rl = 0;
		else if (su==order-1) rl = 1;
		
		switch (rl) {
		case 0:
			int n = freePool.get(su).size();
			if (n==0) return;
			//Check if the residue can be moved, there must exist residues between groups
			if (freePool.get(su).get(n-1) > block.get(su).get(subunitLen-1)){
				freePool.get(su+1).add(0,freePool.get(su).get(n-1));
				freePool.get(su).remove(n-1);
			}
			break;
			
		case 1:
			//Check if the residue can be moved, there must exist residues between groups
			if (freePool.get(su).size() == 0) return;
			if (freePool.get(su).get(0) < block.get(su).get(0)){
				freePool.get(su-1).add(freePool.get(su).get(0));
				freePool.get(su).remove(0);
			}
			break;
		}
	}
	
	/**
	 *  It extends at the beginning or end a group of consecutive aligned residues by moving the residues from the
	 *  freePool to the block.
	 */
	private void expandBlock(){
		
		//Initialize a random generator number
		Random rnd = new Random();
		boolean expandable = true;
		
		//If any freePool is empty, the block cannot be extended (one or more subunits cannot)
		for (int su=0; su<order; su++){
			if (freePool.get(su).size()==0) expandable = false;
		}
		
		if(expandable){
			
			int rl = rnd.nextInt(2);  //Select between expanding right (0) or left (1)
			int res = rnd.nextInt(subunitLen); //Residue as a pivot to expand the subunits
			
			switch (rl) {
			case 0:
				
				//Check that there is at least one residue in the freePool to the right (bigger) than the pivot in each subunit
				for (int su=0; su<order; su++){
					if (freePool.get(su).get(freePool.get(su).size()-1)<block.get(su).get(res)) expandable = false;
				}
				if (!expandable) return;
				
				//Find the next expandable group of residues from the pivot to the right (bigger)
				int rightRes = res+1;
				for (int i=res; i<subunitLen; i++){
					//Break if the end of the aligned residues has been found
					if (rightRes == subunitLen){
						break;
					}
					boolean groupFound = true;
					//Only a group is found when there is a gap in all the subunits in that position
					for (int su=0; su<order; su++){
						if (block.get(su).get(rightRes)-1 == block.get(su).get(i)){
							groupFound = false;
							break;
						}
					}
					if (groupFound){
						break;
					}
					rightRes++;
				}
				
				//Special case: when the rightRes==subunitLen and there is no freePool residue higher, avoid adding a residue outside the subunit
				for (int su=0; su<order; su++){
					if (rightRes==subunitLen && freePool.get(su).get(freePool.get(su).size()-1)<block.get(su).get(subunitLen-1)){
						expandable = false;
					}
				}
				if (!expandable) return;
				
				
				//Expand the block with the residues and delete them from the freePool
				for (int su=0; su<order; su++){
					Integer residue = block.get(su).get(rightRes-1)+1;
					if (rightRes == subunitLen){
						block.get(su).add(residue);
						freePool.get(su).remove(residue);
					} else{
						block.get(su).add(rightRes, residue);
						freePool.get(su).remove(residue);
					}
				}
				subunitLen++;
				break;
				
			case 1:
				
				//Check that there is at least one residue in the freePool to the left (smaller) than the first block in each subunit
				for (int su=0; su<order; su++){
					if (freePool.get(su).get(0)>block.get(su).get(res)) expandable = false;
				}
				if (!expandable) return;
				
				//Find the next expandable group of residues from the pivot to the left (smaller)
				int leftRes = res-1;
				for (int i=res; i>0; i--){
					//Break if the start of the aligned residues has been found
					if (leftRes < 0){
						break;
					}
					boolean groupFound = true;
					//Only a group is found when there is a gap in all the subunits in that position
					for (int su=0; su<order; su++){
						if (block.get(su).get(leftRes)+1 == block.get(su).get(i)){
							groupFound = false;
							break;
						}
					}
					if (groupFound){
						break;
					}
					leftRes--;
				}
				
				//Special case: when the leftRes==-1 and there is no freePool residue lower, avoid adding a residue outside the subunit
				for (int su=0; su<order; su++){
					if (leftRes<0 && freePool.get(su).get(0)>block.get(su).get(0)){
						expandable = false;
					}
				}
				if (!expandable) return;
				
				//Expand the block with the residues and delete them from the freePool
				for (int su=0; su<order; su++){
					Integer residue = block.get(su).get(leftRes+1)-1;
					block.get(su).add(leftRes+1,residue);
					freePool.get(su).remove(residue);
				}
				subunitLen++;
				break;
			}
		}
	}
	
	/**
	 *  It moves aligned residues from the start or end of a group of consecutive aligned residues 
	 *  from the block to the freePool.
	 */
	private void shrinkBlock(){
		
		if (subunitLen < 2) return; //If there is only one aligned residue in the subunits do not shrink.
		
		//Initialize a random generator number
		Random rnd = new Random();
		int rl = rnd.nextInt(2);  //Select between shrink right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to shrink the subunits
		
		switch (rl) {
		case 0:
			//Find the last aligned group of residues before a gap to the right (bigger)
			int rightRes = res+1;
			for (int i=res; i<subunitLen; i++){
				//Break if the end of the aligned residues has been found
				if (rightRes == subunitLen){
					break;
				}
				boolean groupFound = true;
				//Only a group is found when there is a gap in all the subunits in the next position
				for (int su=0; su<order; su++){
					if (block.get(su).get(rightRes)-1 == block.get(su).get(i)){
						groupFound = false;
						break;
					}
				}
				if (groupFound){
					break;
				}
				rightRes++;
			}
			//Shrink the block and add the residues to the freePool
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(rightRes-1);
				block.get(su).remove(rightRes-1);
				freePool.get(su).add(residue);
				Collections.sort(freePool.get(su));
			}
			subunitLen--;
			break;
			
		case 1:
			//Find the first aligned group of residues after a gap to the left (smaller) than the pivot
			int leftRes = res-1;
			for (int i=res; i>0; i--){
				//Break if the start of the aligned residues has been found
				if (leftRes < 0){
					break;
				}
				boolean groupFound = true;
				//Only a group is found when there is a gap in all the subunits in that position
				for (int su=0; su<order; su++){
					if (block.get(su).get(leftRes)+1 == block.get(su).get(i)){
						groupFound = false;
						break;
					}
				}
				if (groupFound){
					break;
				}
				leftRes--;
			}
			
			//Shrink the block and add the residues to the freePool
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(leftRes+1);
				block.get(su).remove(leftRes+1);
				freePool.get(su).add(residue);
				Collections.sort(freePool.get(su));
			}
			subunitLen--;
			break;
		}
	}
	
	private void splitBlock(){
		
		if (subunitLen < 2) return; //If there is only one aligned residue in the subunits do not split and shrink.
		
		//Initialize a random generator number
		Random rnd = new Random();
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to split the subunits
		
		for (int su=0; su<order; su++){
			Integer residue = block.get(su).get(res);
			block.get(su).remove(res);
			freePool.get(su).add(residue);
			Collections.sort(freePool.get(su));
		}
		subunitLen--;
	}
	
	/**
	 *  Calculates the average RMSD of all possible rotation superimpositions of the molecule, corresponding to the
	 *  aligned residues in block. 
	 */
	private void updateAvgRMSDandScore() throws StructureException{
		
		double avgRMSD = 0.0;
		double avgScore = 0.0;
		
		//Construct all possible rotations of the molecule (order-1 possible, index i)
		for (int i=1; i<order; i++){
			RMSDandScore(i);
			avgRMSD += rmsd;
			avgScore += score;
		}
		
		//Construct all possible subunit superimpositions of the molecule: order*(order-1)/2 possible
		/*for (int i=0; i<order; i++){
			for (int j=i+1; j<order; j++){
				subunitRMSDandScore(i, j);
				avgRMSD += rmsd;
				avgScore += score;
			}
		}*/
		
		rmsd = avgRMSD / (order-1);
		score = avgScore / (order-1);
	}
	
	/**
	 *  Raw implementation of the RMSD and Score calculation for a better algorithm efficiency.
	 */
	private void RMSDandScore(int rotation) throws StructureException{
		
		Atom[] arr1 = new Atom[subunitLen*order];
		Atom[] arr2 = new Atom[subunitLen*order];
		int pos = 0;
		
		//Calculate the aligned atom arrays
		for (int j=0; j<order; j++){
			for (int k=0; k<subunitLen; k++){
				arr1[pos] = ca[block.get(j).get(k)];
				arr2[pos] = (Atom) ca[block.get((rotation+j)%order).get(k)].clone();
				pos++;
			}
		}
		
		//Superimpose the two structures in correspondence to the new alignment
		SVDSuperimposer svd = new SVDSuperimposer(arr1, arr2);
		Matrix matrix = svd.getRotation();
		Atom shift = svd.getTranslation();
		
		for (Atom a : arr2) {
			Calc.rotate(a, matrix);
			Calc.shift(a, shift);
		}
		
		//Get the rmsd and score of the rotation
		rmsd = SVDSuperimposer.getRMS(arr1, arr2);
		score = SVDSuperimposer.getTMScore(arr1, arr2, ca.length, ca.length);
	}
	
	/**
	 *  Calculate RMSD and score only considering the subunits to superimpose, not the whole protein.
	 *  NOTE: Did not result to be useful alone, because the optimization converges to parts that are not subunits,
	 *        like big insertions, if they have helices or sheets. Consider it as additional information.
	 */
	private void subunitRMSDandScore(int su1, int su2) throws StructureException{
		
		//Create the atom arrays corresponding to the first and second subunits only
		int su1start = 0;
		int su1end = 0;
		int su2start = 0;
		int su2end = 0;
		
		//su1
		su1start = block.get(su1).get(0);
		su1end = block.get(su1).get(block.get(su1).size()-1)+1;
		
		//su2
		su2start = block.get(su2).get(0);
		su2end = block.get(su2).get(block.get(su2).size()-1)+1;
		/*if (freePool.get(su2).size()==0){
			su2start = block.get(su2).get(0);
			su2end = block.get(su2).get(block.get(su2).size()-1)+1;
		} else {
			//su2 start
			if (block.get(su2).get(0) < freePool.get(su2).get(0)){
				su2start = block.get(su2).get(0);
			} else {
				su2start = freePool.get(su2).get(0);
			}
			//su1 end
			if (block.get(su2).get(block.get(su2).size()-1) > freePool.get(su2).get(freePool.get(su2).size()-1)){
				su2end = block.get(su2).get(block.get(su2).size()-1)+1;
			} else {
				su2end = freePool.get(su2).get(freePool.get(su2).size()-1)+1;
			}
		}*/
		Atom[] ca1block = Arrays.copyOfRange(ca, su1start, su1end);
		Atom[] ca2block = Arrays.copyOfRange(ca, su2start, su2end);
		
		//Extract the aligned atom arrays
		Atom[] aligned1 = new Atom[subunitLen];
		Atom[] aligned2 = new Atom[subunitLen];
		
		for (int k=0; k<subunitLen; k++){
			aligned1[k] = ca[block.get(su1).get(k)];
			aligned2[k] = (Atom) ca[block.get(su2).get(k)].clone();
		}
		
		//Superimpose the two structures in correspondence to the new alignment
		SVDSuperimposer svd = new SVDSuperimposer(aligned1, aligned2);
		Matrix matrix = svd.getRotation();
		Atom shift = svd.getTranslation();
		
		for (Atom a : aligned2) {
			Calc.rotate(a, matrix);
			Calc.shift(a, shift);
		}
		
		//Get the rmsd and score of the rotation
		rmsd = SVDSuperimposer.getRMS(aligned1, aligned2);
		score = SVDSuperimposer.getTMScore(aligned1, aligned2, ca1block.length, ca2block.length);
	}
	
	/**
	 *  Calculate the multiple sequence alignment as strings from the block and freePool residues.
	 */
	private void updateSeqAln() {
	    
		//Store the current positions of the alignment, block and free
		int[] blockPos = new int[order];
		int[] freePos = new int[order];
		
		//Initialize sequence variables
		alnSequences = new String[order];
		Arrays.fill(alnSequences,new String());
		alnSymbols = new String();
		
		while(true){
			boolean stop = true;
			int gaps = 0;
			char[] provisional = new char[order];
		    
			for (int su=0; su<order; su++){
				//If the subunit residues ran out (no more free or block) insert a gap
				if (blockPos[su] == subunitLen && freePos[su] == freePool.get(su).size()){
					provisional[su] = '-';
					gaps++;
				}
				//If there are no more aligned residues increment the gap and add the freePool residues
				else if (blockPos[su] == subunitLen){
					provisional[su] = StructureTools.get1LetterCode(ca[freePool.get(su).get(freePos[su])].getGroup().getPDBName());
					freePos[su]++;
					gaps++;
				}
				//If there are no more free residues insert a provisional dash
				else if (freePos[su] == freePool.get(su).size()){
					provisional[su] = '-';
				}
				//If the next residue is from the freePool, because the residue is lower than the one in the block
				else if (freePool.get(su).get(freePos[su]) < block.get(su).get(blockPos[su])){
					provisional[su] = StructureTools.get1LetterCode(ca[freePool.get(su).get(freePos[su])].getGroup().getPDBName());
					freePos[su]++;
					gaps++;
				}
				//If the next residue is from the block insert a provisional dash
				else {
					provisional[su] = '-';
				}
			}
			//If there are no gaps add the aligned residues
			if (gaps==0){
				alnSymbols += "|";
				for (int su=0; su<order; su++){
					alnSequences[su] += StructureTools.get1LetterCode(ca[block.get(su).get(blockPos[su])].getGroup().getPDBName());
					blockPos[su]++;
				}
			}
			//If any sequence has a gap add the provisional characters
			else {
				alnSymbols += " ";
				for (int su=0; su<order; su++){
					alnSequences[su] += provisional[su];
				}
			}
			
			//Stop if all of the residues in the block and freePool have been analyzed
			for (int q=0; q<order; q++){
				if (freePos[q] != freePool.get(q).size() || blockPos[q] != subunitLen)
					stop = false;
			}
			if (stop) break;
		}
		
		if (debug){
			System.out.println("SEQUENCE ALIGNMENT:");
			for (int su=0; su<order; su++){
				System.out.println(alnSequences[su]);
				if (su!=order-1) System.out.println(alnSymbols);
			}
			System.out.println("Subunit Length: "+subunitLen);
		}
	}
	
	/**
	 * Save the evolution of the optimization process as a csv file.
	 */
	private void saveHistory(String filePath) throws IOException{
		
	    FileWriter writer = new FileWriter(filePath);
	    writer.append("Step,Subunit.Length,RMSD,Score\n");
	    
	    for (int i=0; i<subunitLenHistory.size(); i++){
	    		writer.append(i+","+subunitLenHistory.get(i)+","+rmsdHistory.get(i)+","+scoreHistory.get(i)+"\n");
	    }
	    
	    writer.flush();
	    writer.close();
	}
	
	/**
	 *  Save the multiple sequence alignment result of the optimization as a fasta alignment
	 */
	private void saveSeqAln(String filePath) throws IOException {
		
		FileWriter writer = new FileWriter(filePath);
		
		for (int su=0; su<order; su++){
			writer.append(">Subunit_"+(su+1)+"\n");
			writer.append(alnSequences[su]+"\n");
			//if (su!=order-1) writer.append(alnSymbols+"\n");
		}
		
		writer.flush();
	    writer.close();
	}

	public static void main(String[] args) throws IOException, StructureException{
		
		long startTime = System.currentTimeMillis();
		
		String name = "4i4q";
		AtomCache cache = new AtomCache();
		Atom[] ca1 = cache.getAtoms(name);
		Atom[] ca2 = cache.getAtoms(name);
		
		CeSymm ceSymm = new CeSymm();
		AFPChain afpChain = ceSymm.align(ca1, ca2);
		afpChain.setName1(name);
		
		long alignTime = System.currentTimeMillis() - startTime;
		
		MultipleAFP multipleAFP = new MultipleAFP(afpChain,ca1);
		
		long optimizationTime   = System.currentTimeMillis() - alignTime - startTime;
		
		if (debug){
			System.out.println(">> Alignment time: "+alignTime+" ms");
			System.out.println(">> Optimization time: "+optimizationTime+" ms");
		}
		
	}

	public AFPChain getAfpChain() {
		return afpChain;
	}

	public void setAfpChain(AFPChain afpChain) {
		this.afpChain = afpChain;
	}

	public String[] getAlnSequences() {
		return alnSequences;
	}

	public void setAlnSequences(String[] alnSequences) {
		this.alnSequences = alnSequences;
	}

	public String getAlnSymbols() {
		return alnSymbols;
	}

	public void setAlnSymbols(String alnSymbols) {
		this.alnSymbols = alnSymbols;
	}
}
