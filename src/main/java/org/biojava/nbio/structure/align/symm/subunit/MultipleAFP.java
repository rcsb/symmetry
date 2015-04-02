package org.biojava.nbio.structure.align.symm.subunit;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CeSymm;
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
	
	//List to store the residues aligned, in the block. Dimensions are: [order][residues in block]
	List<ArrayList<Integer>> block;
	//List to store the residues not aligned, that are part of the free pool. Dimensions are: [order][residues in the pool]
	List<ArrayList<Integer>> freePool;
	
	//They store the average RMSD and score of all possible rotation superpositions (alignments) of the protein
	double rmsd;
	double score;
	
	MultipleAFP(AFPChain afp, Atom[] ca1){
		
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
		
		//Add the subunit residues aligned, previously found, to the block
		for (int i=0; i<order; i++){
			ArrayList<Integer> residues = new ArrayList<Integer>();
			for (int j=0; j<afpChain.getOptLen()[i]; j++){
				residues.add(afpChain.getOptAln()[i][0][j]);
				alreadySeen.add(afpChain.getOptAln()[i][0][j]);
			}
			block.add(residues);
			freePool.add(new ArrayList<Integer>());
		}
		int subunitSize = block.get(0).size();
		
		//Add any residue not aligned to the free pool, in the corresponding subunit region (flexible, initial assignment)
		for (int i=0; i<ca.length; i++){
			if (alreadySeen.contains(i)) continue;
			for (int j=0; j<order; j++){
				if (i<block.get(j).get(subunitSize-1) && !alreadySeen.contains(i)){
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
			optimizeMC(10000);
			saveHistory("/scratch/mcopt/1tim.a_MCopt_0.csv");
		} catch (StructureException e) {
			e.printStackTrace();
		} catch (IOException e) {
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
		double C = 1.0;  //Constant for the probability of acceptance for bad moves, will be decreased
		
		int i = 1;
		
		while (i<maxIter && conv<1000){
			
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
			
			if (AS<0){
				
				//Probability of accepting the new alignment
				double prob = Math.min(Math.max((C+10*AS)/(i),0.0),1.0);
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
			
			if (debug) 	System.out.println(i+": --prob: "+Math.min(Math.max((C+10*AS)/(i),0.0),1.0)+", --score: "+AS+", --conv: "+conv);
			
			subunitLenHistory.add(subunitLen);
			rmsdHistory.add(rmsd);
			scoreHistory.add(score);
			
			i++;
		}
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
		boolean expanded = false;
		int iter = 0;
		
		//If any freePool is empty, the block cannot be extended (one or more subunits cannot)
		for (int su=0; su<order; su++){
			if (freePool.get(su).size()==0) expanded = true;
		}
		
		while(!expanded && iter<order){
			
			iter++;
			int rl = rnd.nextInt(2);  //Select between expanding right (0) or left (1)
			int res = rnd.nextInt(subunitLen); //Residue as a pivot to expand the subunits
			
			switch (rl) {
			case 0:
				
				expanded = true;
				//Check that there is at least one residue in the freePool to the right (bigger) than the pivot in each subunit
				for (int su=0; su<order; su++){
					if (freePool.get(su).get(freePool.get(su).size()-1)<block.get(su).get(res)) expanded = false;
				}
				if (!expanded) continue;
				
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
						expanded = false;
					}
				}
				if (!expanded) continue;
				
				
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
				
				expanded = true;
				//Check that there is at least one residue in the freePool to the left (smaller) than the first block in each subunit
				for (int su=0; su<order; su++){
					if (freePool.get(su).get(0)>block.get(su).get(res)) expanded = false;
				}
				if (!expanded) continue;
				
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
						expanded = false;
					}
				}
				if (!expanded) continue;
				
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
	 *  If it is too slow the replaceOptAln method can be changed for a raw implementation of the RMSD calculation.
	 */
	private void updateAvgRMSDandScore() throws StructureException{
		
		double avgRMSD = 0.0;
		double avgScore = 0.0;
		
		//Construct all possible rotations of the molecule (order-1 possible, index i)
		for (int i=1; i<order; i++){
			getRMSDandScore(i);
			avgRMSD += rmsd;
			avgScore += score;
		}
		rmsd = avgRMSD / (order-1);
		score = avgScore / (order-1);
	}
	
	/**
	 *  Raw implementation of the RMSD and Score calculation for a better algorithm efficiency.
	 */
	private void getRMSDandScore(int rotation) throws StructureException{
		
		int n = block.get(0).size();  //length of the subunit
		Atom[] arr1 = new Atom[n*order];
		Atom[] arr2 = new Atom[n*order];
		int pos = 0;
		
		//Calculate the aligned atom arrays
		for (int j=0; j<order; j++){
			for (int k=0; k<n; k++){
				arr1[pos] = ca[block.get(j).get(k)];
				arr2[pos] = (Atom) ca[block.get((rotation+j)%order).get(k)].clone();
				pos++;
			}
		}
		
		//Superimpose the two structures in correspondance to the new alignment
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
	
	private void saveHistory(String filePath) throws IOException{
		
	    FileWriter writer = new FileWriter(filePath);
	    writer.append("Step,Subunit.Length,RMSD,Score\n");
	    
	    for (int i=0; i<subunitLenHistory.size(); i++){
	    		writer.append(i+","+subunitLenHistory.get(i)+","+rmsdHistory.get(i)+","+scoreHistory.get(i)+"\n");
	    }
	    
	    writer.flush();
	    writer.close();
	}
	
	public static void main(String[] args) throws IOException, StructureException{
		
		String name = "1tim.a";
		AtomCache cache = new AtomCache();
		Atom[] ca1 = cache.getAtoms(name);
		Atom[] ca2 = cache.getAtoms(name);
		
		CeSymm ceSymm = new CeSymm();
		AFPChain afpChain = ceSymm.align(ca1, ca2);
		
		MultipleAFP multipleAFP = new MultipleAFP(afpChain,ca1);
		
		System.out.println("Finished!");
		
	}
}
