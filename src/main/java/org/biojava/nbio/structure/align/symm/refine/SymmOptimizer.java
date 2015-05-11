package org.biojava.nbio.structure.align.symm.refine;

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
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Creates a refined alignment by a Monte Carlo optimization of a subunit multiple alignment.
 * 
 * This class pretends to use the CEMC approach for multiple structural alignment of the subunits
 * using the refined pairwise alignment obtained from the align method.
 * 
 * @author lafita
 */
public class SymmOptimizer {

	private static final boolean debug = false;
	
	//Optimization parameters
	private static final int Lmin = 4; //Minimum block length of aligned residues
	public int iterFactor = 100; //Factor to control the max number of iterations of optimization
	public double C = 5; //Probability function constant (probability of acceptance for bad moves)
	
	//Score function parameters
	private static final double M = 20.0; //Maximum score of a match
	private static final double A = 10.0; //Penalty for alignment distances
	public double d0; //Maximum distance that is not penalized - chosen from initial alignment RMSD
	
	private AFPChain afpChain;
	private Atom[] ca;
	private int order;
	public int subunitLen;
	
	//Variables that store the history of the optimization, in order to be able to plot the evolution of the system.
	private List<Integer> subunitLenHistory;
	private List<Double> rmsdHistory;
	private List<Double> scoreHistory;
	
	//List to store the residues aligned, in the block. Dimensions are: [order][subunitLen]
	private List<ArrayList<Integer>> block;
	//List to store the residues not aligned, that are part of the free pool. Dimensions are: [order][residues in the pool]
	private List<ArrayList<Integer>> freePool;
	
	public double rmsd;     // Average RMSD of all rotation superpositions
	public double tmScore;  // Average TM-score of all rotation superpositions
	public double mcScore;  // Optimization score, calculated as the original CEMC algorithm
	double[] colDistances;   //Stores the average distance of the residues in a column. Length: subunitLen
	
	//Multiple alignment String sequences
	String[] alnSequences;
	String alnSymbols;
	
	public SymmOptimizer() {
		super();
	}
	
	public AFPChain optimize(AFPChain seedAFP, Atom[] ca1, Atom[] ca2, int order)
			throws RefinerFailedException,StructureException {
		
		//No multiple alignment can be generated if there is only one subunit.
		if (order == 1) return seedAFP;
		if (afpChain.getBlockNum() == 0) throw new RefinerFailedException("Empty seed alignment");
		
		//Set parameters from initial alignment
		d0 = Math.max(seedAFP.getTotalRmsdOpt()*2,5);
		
		initialize(seedAFP, ca1);
		optimizeMC(iterFactor*ca.length);
		
		/*try {
			saveHistory("/scratch/mcopt/"+afpChain.getName1()+"_MC.csv");
			if (debug) System.out.println("Saved history.");
		} catch (IOException e) {
			e.printStackTrace();
		}*/
		
		return afpChain;
	}
	
	private void initialize(AFPChain afp, Atom[] ca1) throws StructureException{
		
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
		//Set the scores and RMSD of the initial state
		updateScore();
		
		/*//Set the constant d0 to control the bad alignment penalty term (A) from the initial state. Now calculate from originalAFP RMSD.
		calculatePenaltyDistance();
		//Calculate MCscore again with the new d0 parameter
		mcScore = scoreFunctionMC(colDistances);*/
		
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
		
		int i = 1;
		
		while (i<maxIter && conv<(maxIter/10)){
			
			moveResidue();  //At the beginning of each iteration move randomly a residue from the freePool
			
			//Save the state of the system in case the modifications are not favorable
			List<ArrayList<Integer>> lastBlock = new ArrayList<ArrayList<Integer>>();
			List<ArrayList<Integer>> lastFreePool = new ArrayList<ArrayList<Integer>>();
			for (int k=0; k<order; k++){
				lastBlock.add((ArrayList<Integer>) block.get(k).clone());
				lastFreePool.add((ArrayList<Integer>) freePool.get(k).clone());
			}
			double lastScore = mcScore;
			double lastRMSD = rmsd;
			double lastTMscore = tmScore;
			double[] lastColDistances = colDistances.clone();
			
			
			boolean moved = false;
			
			while (!moved){
				//Randomly select one of the steps to modify the alignment
				int move = rnd.nextInt(4);
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
				}
			}
			
			//Get the properties of the new alignment
			updateScore();
			
			double AS = mcScore-lastScore;  //Change in the optimization Score
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
					colDistances = lastColDistances;
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
		
		afpChain = AlignmentTools.replaceOptAln(newAlgn, afpChain, ca, ca);
		updateSeqAln();
	}

	/**
	 *  Move all the block residues of one subunit one position to the left or right and move the corresponding
	 *  boundary residues from the freePool to the block, and viceversa.
	 */
	private boolean shiftRow(){
		
		//Initialize a random generator number
		Random rnd = new Random();
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
		
		//Initialize a random generator number
		Random rnd = new Random();
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
			if ((rightRes-res) <= Lmin) return moved;  //If the block (consecutive residues) is short don't shrink
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
			if ((res-leftRes) <= Lmin) return moved;  //If the block (consecutive residues) is short don't shrink
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
		if (subunitLen <= Lmin) return moved; //Let split moves everywhere if the subunit is larger than the minimum block size
		
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
		moved = true;
		return moved;
	}
	
	/**
	 *  Calculates the average RMSD and Score of all possible rotation superimpositions of the molecule, corresponding to the
	 *  aligned residues in block. It also updates the Monte Carlo score, the optimized value, from the distances matrix.
	 */
	private void updateScore() throws StructureException{
		
		//The first index is the score and the second the RMSD
		double[] ScoreRMSD = {0.0,0.0};
		double[] distances = new double[subunitLen]; //Stores the average distances between the residues in a column (every index corresponds to a column)
		
		//Reset old values
		rmsd = 0.0;
		tmScore = 0.0;
		mcScore = 0.0;
		
		//Construct all possible rotations of the molecule (order-1 possible, index i)
		for (int i=1; i<order; i++){
			calculateScore(i, ScoreRMSD, distances);
			rmsd += ScoreRMSD[1];
			tmScore += ScoreRMSD[0];
		}
		//Divide the colDistances entries for the total number to get the average
		int total = (order)*(order-1);
		for (int i=0; i<subunitLen; i++) distances[i] /= total;
		colDistances = distances;
		
		//Assign the new values to the member variables
		rmsd /= (order-1);
		tmScore /= (order-1);
		mcScore = scoreFunctionMC(colDistances);
	}
	
	/**
	 *  Raw implementation of the RMSD, TMscore and distances calculation for a better algorithm efficiency.
	 */
	private void calculateScore(int rotation, double[] ScoreRMSD, double[] distances) throws StructureException{
		
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
			double distance = Calc.getDistance(arr1[k], arr2[k]);
			distances[k%subunitLen] += distance;
		}
	}
	
	/**
	 *  Calculates the optimization score from the column average distances.
	 *  
	 *  Function: sum(M/(d1/d0)^2)  and add the penalty A if (d1>d0).
	 */
	private double scoreFunctionMC(double[] colDistances){
		
		double score = 0.0;
		
		//Loop through all the columns
		for (int col=0; col<subunitLen; col++){
			
			double d1 = colDistances[col];
			double colScore = M/Math.pow(1+d1/d0,2);
			
			if (d1>d0) colScore-=A;
			score += colScore;
		}
		return score;
	}
	
	/**
	 *  Calculates the probability of accepting a bad move given the iteration step and the score change.
	 *  
	 *  Function: p=(C-AS)/m^0.5   *from the CEMC algorithm.
	 *  Added a normalization factor so that the probability approaches 0 when the maxIter is reached.
	 */
	private double probabilityFunction(double AS, int m, int maxIter) {
		
		double prob = (C+AS)/Math.sqrt(m);
		double norm = (1-(m*1.0)/maxIter);  //Normalization factor
		return Math.min(Math.max(prob*norm,0.0),1.0);
	}
	
	/**
	 *  Calculate the maximum distance which is not penalized in the score function. Only used at the beginning.
	 *  Options are:
	 *    1- Pick the 90% of the distances range (softer condition, results in longer subunits).
	 *    2- Pick the average value of the top 10% distances (can be softer than 1 depending on the range scale).
	 *    3- Pick the value at the boundary of the top 10% distances (hardest condition, restricts the subunits to the core only).
	 *    
	 *  Set a minimum distance to avoid short refined alignments.
	 */
	@SuppressWarnings("unused")
	private void calculatePenaltyDistance(){
	
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
		
		d0=Math.max(d1,4);
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
	    		writer.append(i*100+","+subunitLenHistory.get(i)+","+rmsdHistory.get(i)+","+scoreHistory.get(i)+"\n");
	    }
	    
	    writer.flush();
	    writer.close();
	}
	
	public static void main(String[] args) throws IOException, StructureException, RefinerFailedException{
		
		//Easy cases: 4i4q, 4dou
		//Hard cases: d2vdka_,d1n6dd3, d1n7na1
		//Better MULTIPLE: 2i5i.a
		String[] names = {//"d2vdka_", "d1n6dd3", "d2g02a1", "d1jofd_", "d1kkta_", "d1pbyb_",	"d1ri6a_", "d1tyqc_", "d1xksa_",  //C7
						  //"d1k32f2", "d1okca_", "d1q7fa_", "d1qlga_", "d1uyox_", "d1wp5a_", "d1zxua1", "d2agsa2", "d2ivza1", //C6
						  //"d1ffta_", "d1i5pa2", "d1jlya1", "d1lnsa1", "d1r5za_", "d1ttua3", "d1vmob_", "d1wd3a2", "d2hyrb1", //C3
						  //"d1m1ha1", "d1pexa_", //C4
						  //"d1vkde_", "d2h2na1", "d2jaja_" //C5
						  "d1b0ja2"  //C2
						  };
		for (String name:names){
			
		int order = 2;
			
		System.out.println(name);
		
		AtomCache cache = new AtomCache();
		Atom[] ca1 = cache.getAtoms(name);
		Atom[] ca2 = cache.getAtoms(name);
		
		CeSymm ceSymm = new CeSymm();
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		AFPChain afpChain = ceSymm.align(ca1, ca2);
		
		//Force the order of symmetry that we want
		SymmOptimizer refiner = new SymmOptimizer();
		AFPChain refinedAFP = refiner.optimize(afpChain, ca1, ca2, order);
		
		afpChain.setName1(name);
		afpChain.setName2(name);
		
		SymmetryJmol jmol = new SymmetryJmol(refinedAFP, ca1);
		}
		
		System.out.println("Finished Alaysis!");
	}
}