package org.biojava.nbio.structure.align.symm.refine;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Optimizes an alignment by a Monte Carlo score optimization of the subunit multiple alignment.
 * Implements Callable in order to parallelize multiple optimizations.
 * This class pretends to use the same CEMC approach for multiple structural alignment of the subunits
 * using the refined pairwise alignment.
 * 
 * @author Aleix Lafita
 * 
 */
public class SymmOptimizer implements Callable<AFPChain> {

	private static final boolean debug = true;  //Prints the optimization moves and saves a file with the history in results
	private SymmetryType type;
	private Random rnd;
	
	//Optimization parameters
	private static final int Lmin = 8;   //Minimum subunit length
	private int iterFactor = 100; //Factor to control the max number of iterations of optimization
	private double C = 20; //Probability function constant (probability of acceptance for bad moves)
	
	//Score function parameters
	private static final double M = 20.0; //Maximum score of a match
	private static final double A = 10.0; //Penalty for alignment distances
	private double d0 = 5; //Maximum distance that is not penalized - chosen from seed alignment
	
	//Alignment Information
	private AFPChain afpChain;
	private Atom[] ca;
	private int order;
	public int subunitLen;
	
	//Multiple Alignment Residues
	private List<List<Integer>> block;     //List to store the residues aligned, in the block. Dimensions are: [order][subunitLen]
	private List<List<Integer>> freePool; 	//List to store the residues not aligned. Dimensions are: [order][residues in the pool]
	
	//Score information
	public double rmsd;     // Average RMSD of all rotation superpositions
	public double tmScore;  // Average TM-score of all rotation superpositions
	private double mcScore;  // Optimization score, calculated as the original CEMC algorithm
	
	//Superposition information
	private Matrix rotation;        //The rotation matrix of symmetry
	private Atom translation;       //The translation of the symmetry
	private double[] colDistances;  //Stores the average distance of the algined residues in a column. Length: subunitLen
	private double[] rowDistances;  //Stores the average distance from one subunit to all the others. Similarity measure between subunits.
	
	//Variables that store the history of the optimization, in order to be able to plot the evolution of the system.
	private List<Integer> subunitLenHistory;
	private List<Double> rmsdHistory;
	private List<Double> scoreHistory;
	
	//Multiple alignment String sequences. TODO Consider moving the methods to the new MultipleAlignment DS.
	private String[] alnSequences;
	private String alnSymbols;
	
	/**
	 * Constructor. Initializes all the variables needed for the optimization.
	 * @param seedAFP AFPChain with the symmetry subunits split in blocks.
	 * @param ca1
	 * @param type
	 * @param seed
	 * @throws RefinerFailedException 
	 * @throws StructureException 
	 */
	public SymmOptimizer(AFPChain seedAFP, Atom[] ca1, SymmetryType type, long seed) throws RefinerFailedException, StructureException {
		
		//No multiple alignment can be generated if there is only one subunit
		this.order = seedAFP.getBlockNum();
		if (order == 1) throw new RefinerFailedException("Optimization: Non-Symmetric Seed Alignment.");
		
		this.type = type;
		rnd = new Random(seed);
		
		//Initialize the variables with the seed alignment
		initialize(seedAFP, ca1);
	}
	
	@Override
	public AFPChain call() throws Exception {
		
		optimizeMC(iterFactor*ca.length);
		
		//Save the history to the results folder in the symmetry project
		try {
			if (debug) saveHistory("src/main/java/results/SymmOptimizerHistory.csv");
		} catch(FileNotFoundException e) {}
		return afpChain;
	}
	
	private void initialize(AFPChain afp, Atom[] ca1) throws StructureException {
		
		//Initialize member variables
		afpChain = afp;
		ca = ca1;
		order = afpChain.getBlockNum();
		subunitLen = afpChain.getOptLen()[0];
		rotation = afpChain.getBlockRotationMatrix()[0];
		translation = afpChain.getBlockShiftVector()[0];
		
		//Initialize alignment variables
		block = new ArrayList<List<Integer>>();
		freePool = new ArrayList<List<Integer>>();
		
		//Store the residues that have been added to the block
		List<Integer> alreadySeen = new ArrayList<Integer>();
		
		//Generate the initial state of the system from the aligned blocks of the AFPChain
		for (int i=0; i<order; i++){
			List<Integer> residues = new ArrayList<Integer>();
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
		//Move randomly residues of the free Pool to ensure that they are distributed equally to all subunits
		for (int i=0; i<ca.length-(subunitLen*order); i++) moveResidue();
		
		
		//Set the scores and RMSD of the initial state (seed alignment). Check the order and 
		switch (type){
		case CLOSED: 
			updateClosedScore();
			calculatePenaltyDistance();
			updateClosedScore();
			break;
		case OPEN:
			//checkOrder();
			updateOpenScore();
			calculatePenaltyDistance();
			updateOpenScore();
			break;
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
	private void optimizeMC(int maxIter) throws StructureException{
		
		//Initialize the history variables
		subunitLenHistory = new ArrayList<Integer>();
		rmsdHistory = new ArrayList<Double>();
		scoreHistory = new ArrayList<Double>();
		
		int conv = 0;  //Number of steps without an alignment improvement
		int i = 1;
		
		while (i<maxIter && conv<(maxIter/20)){
			
			moveResidue();  //At the beginning of each iteration move randomly a residue from the freePool
			
			//Save the state of the system in case the modifications are not favorable
			List<List<Integer>> lastBlock = new ArrayList<List<Integer>>();
			List<List<Integer>> lastFreePool = new ArrayList<List<Integer>>();
			for (int k=0; k<order; k++){
				List<Integer> b = new ArrayList<Integer>();
				List<Integer> p = new ArrayList<Integer>();
				for (int l=0; l<subunitLen; l++) b.add(block.get(k).get(l));
				for (int l=0; l<freePool.get(k).size(); l++) p.add(freePool.get(k).get(l));
				lastBlock.add(b);
				lastFreePool.add(p);
			}
			double lastScore = mcScore;
			double lastRMSD = rmsd;
			double lastTMscore = tmScore;
			
			boolean moved = false;
			
			while (!moved){
				//Randomly select one of the steps to modify the alignment
				int move = rnd.nextInt(5);
				switch (move){
				case 0: moved = shiftRow();
						if (debug) System.out.println("did shift");
						break;
				case 1: moved = expandBlock();
						if (debug) System.out.println("did expand");
						break;
				case 2: moved = shrinkBlock();
						if (debug) System.out.println("did shrink");
						break;
				case 3: moved = splitBlock();
						if (debug) System.out.println("did split");
						break;
				case 4: moved = shiftAll();
						if (debug) System.out.println("did shift all");
						break;
				}
			}
			
			//Get the properties of the new alignment
			switch (type){
			case CLOSED: 
				updateClosedScore();
				break;
			case OPEN: 
				updateOpenScore();
				break;
			}
			
			//Calculate change in the optimization Score
			double AS = mcScore-lastScore;
			double prob=1.0;
			
			if (AS<0){
				
				//Probability of accepting the new alignment given that produces a negative score change
				prob = probabilityFunction(AS,i,maxIter);
				double p = rnd.nextDouble();
				//Reject the move
				if (p>prob){
					block = lastBlock;
					freePool = lastFreePool;
					subunitLen = block.get(0).size();
					mcScore = lastScore;
					rmsd = lastRMSD;
					tmScore = lastTMscore;
					conv ++; //Increment the number of steps without a change in score
					
				} else conv = 0;
				
			} else conv=0;
			
			if (debug) 	System.out.println(i+": --prob: "+prob+", --score: "+AS+", --conv: "+conv);
			
			if (i%100==1){
				subunitLenHistory.add(subunitLen);
				rmsdHistory.add(rmsd);
				scoreHistory.add(mcScore);
			}
			
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
		
		//Generate the optimized AFPChain to return - override superimposition information
		afpChain = AlignmentTools.replaceOptAln(newAlgn, afpChain, ca, ca);
		Arrays.fill(afpChain.getBlockRotationMatrix(),rotation);
		Arrays.fill(afpChain.getBlockShiftVector(),translation);
		afpChain.setTMScore(tmScore);
		afpChain.setTotalRmsdOpt(rmsd);
		afpChain.setAlignScore(mcScore);
		
		updateSeqAln();
	}

	/**
	 *  Move all the block residues of one subunit one position to the left or right and move the corresponding
	 *  boundary residues from the freePool to the block, and viceversa.
	 */
	private boolean shiftRow(){
		
		boolean moved = false;

		int su = rnd.nextInt(order); //Select randomly the subunit that is going to be shifted
		int rl = rnd.nextInt(2);  //Select between moving right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to make the shift
		
		if (freePool.get(su).size()==0) return moved;  //If the freePool is empty the subunit cannot be shifted
		
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
				
				moved = true;
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
				
				moved = true;
			}
			break;
		}
		return moved;
	}
	
	/**
	 *  Move all the subunits one position to the left or right. Corresponds to shrinking all subunits at
	 *  one extreme and extending one position at the other. This move is specially designed for symmetry
	 *  because the boundaries of the subunits are not defined and they need to be moved as well.
	 */
	private boolean shiftAll(){
		
		boolean moved = false;
		int rl = rnd.nextInt(2);  //Select between moving right (0) or left (1)
			
		switch(rl){
		case 0: //Move subunits to the right
			
			//There needs to be at least one residue available to the right (end) to make the shift
			if (freePool.get(order-1).size() == 0) return moved;  //The freePool is empty
			else if (freePool.get(order-1).get(freePool.get(order-1).size()-1) < block.get(order-1).get(subunitLen-1)){
				//If it is true it means that there are no available residues at the rightmost of the last subunit
				return moved;
			}
			//This list will store all the residues added to the subunits (needed to maintain the freePool correct)
			List<Integer> extendedResiduesR = new ArrayList<Integer>();
			
			//Loop through all the subunits from last to first
			for (int su=order-1; su>=0; su--){
				//Extend the subunit to the right
				Integer extendRes = block.get(su).get(subunitLen-1)+1;
				extendedResiduesR.add(extendRes);
				Integer shrinkRes = block.get(su).get(0);
				block.get(su).add(extendRes);
				freePool.get(su).add(shrinkRes);
				block.get(su).remove(shrinkRes);
				Collections.sort(freePool.get(su));
			}
			for (Integer res:extendedResiduesR){
				//Remove the residues used for extension from the freePool
				for (int su=0; su<order; su++){
					if (freePool.get(su).remove(res)) break;
				}
			}
			moved = true;
			break;
			
		case 1: //Move subunits to the left
			
			//There needs to be at least one residue available to the left (start) to make the shift
			if (freePool.get(0).size() == 0) return moved;  //The freePool is empty
			else if (freePool.get(0).get(0) > block.get(0).get(0)){
				//If it is true it means that there are no available residues at the start of the first subunit
				return moved;
			}
			//This list will store all the residues added to the subunits (needed to maintain the freePool correct)
			List<Integer> extendedResiduesL = new ArrayList<Integer>();
			
			//Loop through all the subunits from last to first
			for (int su=order-1; su>=0; su--){
				//Extend the subunit to the left
				Integer extendRes = block.get(su).get(0)-1;
				extendedResiduesL.add(extendRes);
				Integer shrinkRes = block.get(su).get(subunitLen-1);
				block.get(su).add(0,extendRes);
				freePool.get(su).add(shrinkRes);
				block.get(su).remove(shrinkRes);
				Collections.sort(freePool.get(su));
			}
			for (Integer res:extendedResiduesL){
				//Remove the residues used for extension from the freePool
				for (int su=0; su<order; su++){
					if (freePool.get(su).remove(res)) break;
				}
			}
			moved = true;
			break;
		}
		return moved;
	}
	
	/**
	 * Move a residue of the freePool from the start of one subunit to the end of the previous, or
	 *                                from the end of one subunit to the start of the next.
	 * Exclude the first and last residues, because makes no sense to move them.
	 * This move does not cause a score change, so it is made each iteration with combination to other moves.
	 */
	private void moveResidue(){
		
		int su = rnd.nextInt(order);
		int rl = rnd.nextInt(2);  //Randomly choose to move one residue from the right-end (0) or left-start (1)
		
		//Exclude first and last residues from being moved.
		if (su==0) rl = 0;
		else if (su==order-1) rl = 1;
		
		switch (rl) {
		case 0:
			int n = freePool.get(su).size();
			if (n==0) return;  //Don't do anything if empty freePool
			//Check if the residue can be moved, there must exist residues between groups
			if (freePool.get(su).get(n-1) > block.get(su).get(subunitLen-1)){
				freePool.get(su+1).add(0,freePool.get(su).get(n-1));
				freePool.get(su).remove(n-1);
			}
			break;
			
		case 1:
			//Check if the residue can be moved, there must exist residues between groups
			if (freePool.get(su).size() == 0) return;  //Don't do anything if empty freePool
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
	private boolean expandBlock(){
		
		boolean moved = false;
		
		//If any freePool is empty, the block cannot be expanded (one or more subunits cannot)
		for (int su=0; su<order; su++){
			if (freePool.get(su).size()==0) return moved;
		}
			
		int rl = rnd.nextInt(2);  //Select between expanding right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to expand the subunits
		
		switch (rl) {
		case 0:
			
			//Check that there is at least one residue in the freePool to the right (bigger) than the pivot in each subunit
			for (int su=0; su<order; su++){
				if (freePool.get(su).get(freePool.get(su).size()-1)<block.get(su).get(res)) return moved;
			}
			
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
					return moved;
				}
			}			
			
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
			moved = true;
			break;
			
		case 1:
			
			//Check that there is at least one residue in the freePool to the left (smaller) than the first block in each subunit
			for (int su=0; su<order; su++){
				if (freePool.get(su).get(0)>block.get(su).get(res)) return moved;
			}
			
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
					return moved;
				}
			}
			
			//Expand the block with the residues and delete them from the freePool
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(leftRes+1)-1;
				block.get(su).add(leftRes+1,residue);
				freePool.get(su).remove(residue);
			}
			subunitLen++;
			moved = true;
			break;
		}
		return moved;
	}
	
	/**
	 *  It moves aligned residues from the start or end of a group of consecutive aligned residues 
	 *  from the block to the freePool.
	 */
	private boolean shrinkBlock(){
		
		boolean moved = false;
		if (subunitLen <= Lmin) return moved; //Do not let shrink moves if the subunit is smaller than the minimum length
		
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
			//if ((rightRes-res) <= AFPmin) return moved;  //If the block (consecutive residues) is short don't shrink
			//Shrink the block and add the residues to the freePool
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(rightRes-1);
				block.get(su).remove(rightRes-1);
				freePool.get(su).add(residue);
				Collections.sort(freePool.get(su));
			}
			subunitLen--;
			moved = true;
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
			//if ((res-leftRes) <= AFPmin) return moved;  //If the block (consecutive residues) is short don't shrink
			//Shrink the block and add the residues to the freePool
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(leftRes+1);
				block.get(su).remove(leftRes+1);
				freePool.get(su).add(residue);
				Collections.sort(freePool.get(su));
			}
			subunitLen--;
			moved = true;
			break;
		}
		return moved;
	}
	
	private boolean splitBlock(){
		
		boolean moved = false;
		if (subunitLen <= Lmin) return moved; //Let split moves everywhere if the subunit is larger than the minimum length
		
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to split the subunits
		
		for (int su=0; su<order; su++){
			Integer residue = block.get(su).get(res);
			block.get(su).remove(res);
			freePool.get(su).add(residue);
			Collections.sort(freePool.get(su));
		}
		subunitLen--;
		moved = true;
		return moved;
	}
	
	/**
	 *  Calculates the average RMSD and Score of all subunit superimpositions of the structure, corresponding to the
	 *  aligned residues in block. It also updates the Monte Carlo score, the optimized value, from the distances matrix.
	 *  It uses a different transformation for every subunit pair (flexible).
	 */
	@Deprecated
	private void updateFlexibleScore() throws StructureException{
		
		//The first index is the score and the second the RMSD
		double[] ScoreRMSD = {0.0,0.0};
		colDistances = new double[subunitLen];
		rowDistances = new double[order];
		
		//Reset old values
		rmsd = 0.0;
		tmScore = 0.0;
		mcScore = 0.0;
		
		//Construct all possible rotations of the molecule (order-1 possible, index i)
		for (int i=0; i<order; i++){
			for (int j=0; j<order; j++){
				if (i!=j){
					calculateFlexibleScore(i, j, ScoreRMSD);
					rmsd += ScoreRMSD[1];
					tmScore += ScoreRMSD[0];
				}
			}
		}
		//Divide the colDistances entries for the total number of comparisons to get the average
		int total = (order)*(order-1);
		for (int i=0; i<subunitLen; i++) colDistances[i] /= total;
		for (int i=0; i<order; i++) rowDistances[i] /= (order-1)*subunitLen;
		
		//Assign the new values to the member variables
		rmsd /= total;
		tmScore /= total;
		mcScore = scoreFunctionCEMC();
		tmScore = scoreFunctionSymmetry();
	}
	
	/**
	 *  Calculates the average RMSD and Score of all subunit superimpositions of the structure, corresponding to the
	 *  aligned residues in block. It also updates the Monte Carlo score, the optimized value, from the distances.
	 *  It uses the same transformation for every subunit pair (obtained from the 2-end vs 1-(end-1) superposition, open symm).
	 */
	private void updateOpenScore() throws StructureException {
		
		//Reset values
		rmsd = 0.0;
		tmScore = 0.0;
		mcScore = 0.0;
		colDistances = new double[subunitLen];
		//rowDistances = new double[order];
		
		//Calculate the aligned atom arrays
		Atom[] arr1 = new Atom[subunitLen*(order-1)];
		Atom[] arr2 = new Atom[subunitLen*(order-1)];
		int pos = 0;
		for (int j=0; j<(order-1); j++){
			for (int k=0; k<subunitLen; k++){
				arr1[pos] = ca[block.get(j).get(k)];
				arr2[pos] = (Atom) ca[block.get(j+1).get(k)].clone();
				pos++;
			}
		}
		//Calculate the new transformation information
		SVDSuperimposer svd = new SVDSuperimposer(arr1, arr2);
		rotation = svd.getRotation();
		translation = svd.getTranslation();
		
		//Construct all possible transformations of the molecule (order-1) and compare pairwise subunits
		for (int i=0; i<order-1; i++){
			//Rotate one more time the atoms of the first alignment array
			Calc.rotate(arr2, rotation);
			Calc.shift(arr2, translation);
			//if (i==0) tmScore += SVDSuperimposer.getTMScore(arr1, arr2, ca.length, ca.length);
			
			for (int j=0; j<order-1-i; j++){
				//Calculate the subunit alignment pairs
				Atom[] su1 = Arrays.copyOfRange(arr1,(j)*subunitLen,(j+1)*subunitLen);
				Atom[] su2 = Arrays.copyOfRange(arr2,(j+i)*subunitLen,(j+i+1)*subunitLen);
				
				//Update score information and distances
				rmsd += SVDSuperimposer.getRMS(su1, su2);
				for (int k=0; k<subunitLen; k++) {
					double distance = Math.abs(Calc.getDistance(su1[k], su2[k]));
					colDistances[k] += distance;
					//rowDistances[j] += distance;
					//rowDistances[(j+i+1)] += distance;
				}
			}
		}
		//Divide the variables for the total number of comparisons to get the average
		int total = ((order)*(order-1))/2;
		for (int i=0; i<subunitLen; i++) colDistances[i] /= total;
		//for (int i=0; i<order; i++) rowDistances[i] /= (order-1)*subunitLen;
		rmsd /= total;
		mcScore = scoreFunctionCEMC();
		tmScore = scoreFunctionSymmetry();
	}
	
	/**
	 *  Calculates the average RMSD and Score of all subunit superimpositions of the structure rotations, corresponding to the
	 *  aligned residues in block. It also updates the Monte Carlo score, the optimized value, from the distances.
	 *  It uses the same transformation for every subunit pair (obtained from the one rotation superposition, closed).
	 *  
	 *  The only difference between the closed and non-closed optimization is that the superposition is made with the whole
	 *  rotation in the closed symmetry, and without the first and last subunits in the non-closed, because they do not align.
	 */
	private void updateClosedScore() throws StructureException{
		
		//Reset values
		rmsd = 0.0;
		tmScore = 0.0;
		mcScore = 0.0;
		colDistances = new double[subunitLen];
		//rowDistances = new double[order];
		
		//Calculate the aligned atom arrays
		Atom[] arr1 = new Atom[subunitLen*order];
		Atom[] arr2 = new Atom[subunitLen*order];
		int pos = 0;
		for (int j=0; j<order; j++){
			for (int k=0; k<subunitLen; k++){
				arr1[pos] = ca[block.get(j).get(k)];
				arr2[pos] = (Atom) ca[block.get((j+1)%order).get(k)].clone();
				pos++;
			}
		}
		//Calculate the new transformation information
		SVDSuperimposer svd = new SVDSuperimposer(arr1, arr2);
		rotation = svd.getRotation();
		translation = svd.getTranslation();
		
		//Construct all possible transformations of the molecule (order-1) and compare pairwise subunits
		for (int i=0; i<order-1; i++){
			//Rotate one more time the atoms of the first atom array
			Calc.rotate(arr2, rotation);
			Calc.shift(arr2, translation);
			//if (i==0) tmScore += SVDSuperimposer.getTMScore(arr1, arr2, ca.length, ca.length);
			
			for (int j=0; j<order-1-i; j++){
				//Calculate the subunit alignment pairs
				Atom[] su1 = Arrays.copyOfRange(arr1,(j)*subunitLen,(j+1)*subunitLen);
				Atom[] su2 = Arrays.copyOfRange(arr2,(j+i)*subunitLen,(j+i+1)*subunitLen);
				
				//Generate the distance information
				rmsd += SVDSuperimposer.getRMS(su1, su2);
				for (int k=0; k<subunitLen; k++) {
					double distance = Math.abs(Calc.getDistance(su1[k], su2[k]));
					colDistances[k] += distance;  //provisional max distances
					//rowDistances[j] += distance;
					//rowDistances[(j+i+1)%order] += distance;
				}
			}
		}
		//Divide the variables for the total number of comparisons to get the average
		int total = ((order)*(order-1))/2;
		for (int i=0; i<subunitLen; i++) colDistances[i] /= total;
		//for (int i=0; i<order; i++) rowDistances[i] /= (order-1)*subunitLen;
		rmsd /= total;
		mcScore = scoreFunctionCEMC();
		tmScore = scoreFunctionSymmetry();
	}
	
	/**
	 *  Calculates the average RMSD and Score of all possible rotation superimpositions of the molecule, corresponding to the
	 *  aligned residues in block. It also updates the Monte Carlo score, the optimized value, from the distances matrix.
	 */
	@Deprecated
	private void updateClosedScoreOld() throws StructureException{
		
		//The first index is the score and the second the RMSD
		double[] ScoreRMSD = {0.0,0.0};
		colDistances = new double[subunitLen]; //Stores the average distances between the residues in a column (every index corresponds to a column)
		
		//Reset old values
		rmsd = 0.0;
		tmScore = 0.0;
		mcScore = 0.0;
		
		//Construct all possible rotations of the molecule (order-1 possible, index i)
		for (int i=1; i<order; i++){
			calculateClosedScore(i, ScoreRMSD);
			rmsd += ScoreRMSD[1];
			tmScore += ScoreRMSD[0];
		}
		//Divide the colDistances entries for the total number to get the average
		int total = (order)*(order-1);
		for (int i=0; i<subunitLen; i++) colDistances[i] /= total;
		for (int i=0; i<order; i++) rowDistances[i] /= (order-1)*subunitLen;
		
		//Assign the new values to the member variables
		rmsd /= (order-1);
		tmScore /= (order-1);
		mcScore = scoreFunctionCEMC();
	}
	
	/**
	 *  Raw implementation of the RMSD, TMscore and distances calculation for a better algorithm efficiency.
	 *  Superimpose flexibly two subunits.
	 */
	@Deprecated
	private void calculateFlexibleScore(int su1, int su2, double[] ScoreRMSD) throws StructureException{
		
		Atom[] arr1 = new Atom[subunitLen];
		Atom[] arr2 = new Atom[subunitLen];
		int pos = 0;
		
		//Calculate the aligned atom arrays
		for (int k=0; k<subunitLen; k++){
			arr1[pos] = ca[block.get(su1).get(k)];
			arr2[pos] = (Atom) ca[block.get(su2).get(k)].clone();
			pos++;
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
		ScoreRMSD[1] = SVDSuperimposer.getRMS(arr1, arr2);
		ScoreRMSD[0] = SVDSuperimposer.getTMScore(arr1, arr2, ca.length, ca.length);
		
		//Calculate the distances between C alpha atoms of the same column and store them in colDistances
		for (int k=0; k<arr1.length; k++){
			double distance = Math.abs(Calc.getDistance(arr1[k], arr2[k]));
			colDistances[k] += distance;
			rowDistances[su1] += distance;
			rowDistances[su2] += distance;
		}
	}
	
	/**
	 *  Raw implementation of the RMSD, TMscore and distances calculation for a better algorithm efficiency.
	 *  Superimpose the original structure with a specified rotation of itself.
	 */
	@Deprecated
	private void calculateClosedScore(int rotation, double[] ScoreRMSD) throws StructureException{
		
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
		ScoreRMSD[1] = SVDSuperimposer.getRMS(arr1, arr2);
		ScoreRMSD[0] = SVDSuperimposer.getTMScore(arr1, arr2, ca.length, ca.length);
		
		//Calculate the distances between C alpha atoms of the same column and store them in colDistances
		for (int k=0; k<arr1.length; k++){
			double distance = Math.abs(Calc.getDistance(arr1[k], arr2[k]));
			colDistances[k%subunitLen] += distance;
			rowDistances[k%order] += distance;
			rowDistances[(k+rotation)%order] += distance;
		}
	}
	
	/**
	 * Checks if there is any subunit/repeat far away from the others and deletes it, decreasing the order.
	 * @throws StructureException 
	 */
	private void checkOrder() throws StructureException{
		
		//First update the scores again because the number of subunits might change distances
		updateOpenScore();
		
		int index = 0;
		double avgDist = 0.0;
		
		//Calculate the index of the subunit with the maximum distance to the others
		for (int i=0; i<order; i++){
			avgDist += rowDistances[i];
			if (rowDistances[index] < rowDistances[i]) index = i;
		}
		avgDist -= rowDistances[index];
		avgDist /= order-1;
		
		//Delete the least similar subunit if its distance to the others is 1.5 times the average (and not below 2.25 A)
		if (rowDistances[index] > 1.5*avgDist && rowDistances[index] > 2.25){
			freePool.get(index).addAll(block.get(index));
			if (index == 0) {
				freePool.get(index+1).addAll(freePool.get(index));
				Collections.sort(freePool.get(index+1));
			} else {
				freePool.get(index-1).addAll(freePool.get(index));
				Collections.sort(freePool.get(index-1));
			}
			freePool.remove(index);
			block.remove(index);
			order -= 1;
			checkOrder(); //call recursively the method to delete all non-relevant subunits (2 is the minimum)
		}
	}
	
	/**
	 *  Calculates the optimization score from the column average distances.
	 *  Function: sum(M/(d1/d0)^2)  and add the penalty A if (d1>d0).
	 */
	private double scoreFunctionCEMC(){
		
		double score = 0.0;
		
		//Loop through all the columns
		for (int col=0; col<subunitLen; col++){
			
			double d1 = colDistances[col];
			
			//CEMC condition
			double colScore = M/(1+(d1*d1)/(d0*d0));
			//if (d1>d0) colScore-=A;
			
			//P2 condition
			colScore -= A;
			colScore *= 2;

			score += colScore;
		}
		return score;
	}
	
	/**
	 *  Calculates a normalized score for the symmetry, to be able to compare between results and set a threshold.
	 *  Function: (1/Ln)*sum(1/(d1/d0)^2). Ln is the max possible length of a subunit (protein length over order).
	 *  The difference with TM score is the normalization factor Ln and the use of average distances.
	 */
	private double scoreFunctionSymmetry(){
		
		double score = 0.0;
		//d0 is calculated as in the TM-score
		double d0 =  1.24 * Math.cbrt((ca.length) - 15.) - 1.8;
		
		//Loop through all the columns
		for (int col=0; col<subunitLen; col++){
			double d1 = colDistances[col];
			double colScore = 1/(1+((d1*d1)/(d0*d0)));
			score += colScore;
		}
		int Ln = ca.length/order;
		
		return score/Ln;
	}
	
	/**
	 *  Calculates the probability of accepting a bad move given the iteration step and the score change.
	 *  
	 *  Function: p=(C-AS)/m   *from the CEMC algorithm.
	 *  Added a normalization factor so that the probability approaches 0 when the maxIter is reached.
	 */
	private double probabilityFunction(double AS, int m, int maxIter) {
		
		double prob = (C+AS)/(m);
		double norm = (1-(m*1.0)/maxIter);  //Normalization factor (step/maxIter)
		return Math.min(Math.max(prob*norm,0.0),0.2);
	}
	
	/**
	 *  Calculate the maximum distance which is not penalized in the score function. Only used at the beginning.
	 *  Options are:
	 *    1- Pick the 90% of the distances range (softer condition, results in longer subunits).
	 *    2- Pick the average value of the top 10% distances (can be softer than 1 depending on the range scale).
	 *    3- Pick the value at the boundary of the top 10% distances (hardest condition, restricts the subunits to the core only).
	 *    4- A function of the RMSD of the seed alignment and the order (this is the softest condition of all, but longer subunits are obtained).
	 *    		Justification: the more subunits the higher the variability between the columns.
	 *    5- A function of the length of the protein (TM-score function)
	 *    
	 *  A minimum distance of 5A is set always to avoid short alignments in the very good symmetric cases.
	 */
	private void calculatePenaltyDistance() {
		
		double[] distances = colDistances.clone();
		Arrays.sort(distances);
		
		//Option 1: 90% of distances range
		double range = distances[distances.length-1] - distances[0];
		double d1 = distances[0] + range*0.9;
		
		//Option 2: average of the top 10% distances
		int index10 = (distances.length-1) - distances.length/10;
		double d2 = 0.0;
		for (double dist:Arrays.copyOfRange(distances,index10,distances.length)) d2+=dist;
		d2 /= (distances.length-index10);
		
		//Option 3: boundary of top 10%
		double d3 = distances[index10];
		
		//Option 4: a function of the seed RMSD
		double d4 = rmsd*(0.5*order);
		
		//Option 5: TM-score function, depends on the length of the protein only
		double d5 =  1.24 * Math.cbrt((ca.length) - 15.) - 1.8;
		
		//d0=d4;
	}
	
	/**
	 *  Calculate the multiple sequence alignment as strings from the block and freePool residues.
	 */
	private void updateSeqAln() {
	    
		//Store the current positions of the alignment, block and freePool
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
	private void saveHistory(String filePath) throws IOException {
		
	    FileWriter writer = new FileWriter(filePath);
	    writer.append("Step,Length,RMSD,Score\n");
	    
	    for (int i=0; i<subunitLenHistory.size(); i++){
	    		writer.append(i*100+","+subunitLenHistory.get(i)+","+rmsdHistory.get(i)+","+scoreHistory.get(i)+"\n");
	    }
	    
	    writer.flush();
	    writer.close();
	}
	
	public static void main(String[] args) throws Exception{
		
		//Easy cases: 4i4q, 4dou
		//Hard cases: d2vdka_,d1n6dd3, d1n7na1
		//Better MULTIPLE: 2i5i.a
		String[] names = {//"d2vdka_", "d1n6dd3", "d2g02a1", "d1jofd_", "d1kkta_", "d1pbyb_",	"d1ri6a_", "d1tyqc_", "d1xksa_",  //C7
						  //"d1k32f2", "d1okca_", "d1q7fa_", "d1qlga_", "d1uyox_", "d1wp5a_", "d1zxua1", "d2agsa2", "d2ivza1", //C6
						  //"d1ffta_", "d1i5pa2", "d1jlya1", "d1lnsa1", "d1r5za_", "d1ttua3", "d1vmob_", "d1wd3a2", "d2hyrb1", //C3
						  //"d1m1ha1", "d1pexa_", //C4
						  //"d1vkde_", "d2h2na1", "d2jaja_", //C5
						  //"d3d5rc1", "d1fz9f_", "d1u4qb3", "d1viia_", "d2a1bf_", "d2axth1",  //C1 border cases
						  //"d1g73b_", "d3pmra_", "d3b5zd2", "d1o5hb_", "d2c2lb1", "d2g38b1", "d1s2xa_", "d1r8ia_"     //SLIP alignment
						  //"d1osya_", "d1xq4a_", "d1yioa2", "d1dcea2", //C2
						  //coiled coil: LPP-56 1jcd, CspB 
						  "4gcr"};  //other
		
		for (String name:names){
			
			System.out.println(name);
			
			AtomCache cache = new AtomCache();
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			CeSymm ceSymm = new CeSymm();
			CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
			params.setSymmetryType(SymmetryType.AUTO);
			params.setRefineMethod(RefineMethod.SINGLE);
			params.setOptimization(true);
			params.setSeed((int) System.currentTimeMillis());
			
			AFPChain afpChain = ceSymm.align(ca1, ca2);
			afpChain.setName1(name);
			afpChain.setName2(name);
			
			SymmetryJmol jmol = new SymmetryJmol(afpChain, ca1);
			//StructureAlignmentJmol jmol2 = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			jmol.setTitle(name);
		}
		
		System.out.println("Finished Alaysis!");
	}
}