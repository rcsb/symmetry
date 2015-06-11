package org.biojava.nbio.structure.align.symm.refine;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.StructureAlignmentException;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Optimizes an alignment by a Monte Carlo score optimization of the subunit multiple alignment.
 * Implements Callable in order to parallelize multiple optimizations.
 * This class pretends to use the same CEMC approach for multiple structural alignment of the subunits
 * using the refined pairwise alignment.
 * Becuase gaps are allowed in some of the subunits, a {@link MultipleAlignment} format has to be returned.
 * 
 * @author Aleix Lafita
 * 
 */
public class SymmGapOptimizer implements Callable<MultipleAlignment> {

	private static final boolean debug = true;  //Prints the optimization moves and saves a file with the history in results
	private SymmetryType type;
	private Random rnd;
	
	//Optimization parameters
	private static final int Rmin = 2;   //Minimum number of aligned structures without a gap
	private static final int Lmin = 8;   //Minimum subunit length
	private int iterFactor = 1000; //Factor to control the max number of iterations of optimization
	private double C = 20; //Probability function constant (probability of acceptance for bad moves)
	
	//Score function parameters
	private static final double M = 20.0; //Maximum score of a match
	private static final double A = 10.0; //Penalty for alignment distances
	private static final double G = 15.0; //Penalty for gaps
	private double d0 = 5; //Maximum distance that is not penalized - chosen from seed alignment
	
	//Alignment Information
	private MultipleAlignment msa;
	private AFPChain seedAFP;
	private Atom[] ca;
	private int order;
	public int subunitLen;
	private int gaps;    //the number of gaps in the alignment
	
	//Multiple Alignment Residues
	private List<List<Integer>> block;     //List to store the residues aligned, in the block. Dimensions are: [order][subunitLen]
	private List<List<Integer>> freePool; 	//List to store the residues not aligned. Dimensions are: [order][residues in the pool]
	
	//Score information
	public double rmsd;     // Average RMSD of all rotation superpositions
	public double tmScore;  // Average TM-score of all rotation superpositions
	private double mcScore;  // Optimization score, calculated as the original CEMC algorithm
	
	//Superposition information
	private Matrix4d transformation;        //The transformation 4D matrix of symmetry
	private double[] colDistances;  //Stores the average distance of the algined residues in a column. Length: subunitLen
	
	//Variables that store the history of the optimization, in order to be able to plot the evolution of the system.
	private List<Integer> subunitLenHistory;
	private List<Double> rmsdHistory;
	private List<Double> scoreHistory;
	
	/**
	 * Constructor. Initializes all the variables needed for the optimization.
	 * @param seedAFP AFPChain with the symmetry subunits split in blocks.
	 * @param ca1
	 * @param type
	 * @param seed
	 * @throws RefinerFailedException 
	 * @throws StructureException 
	 * @throws StructureAlignmentException 
	 */
	public SymmGapOptimizer(AFPChain seedAFP, Atom[] ca1, SymmetryType type, long seed) throws RefinerFailedException, StructureException, StructureAlignmentException {
		
		//No multiple alignment can be generated if there is only one subunit
		this.order = seedAFP.getBlockNum();
		if (order == 1) throw new RefinerFailedException("Optimization: Non-Symmetric Seed Alignment.");
		
		this.type = type;
		rnd = new Random(seed);
		
		//Initialize the variables with the seed alignment
		initialize(seedAFP, ca1);
	}
	
	@Override
	public MultipleAlignment call() throws Exception {
		
		optimizeMC(iterFactor*ca.length);
		
		//Save the history to the results folder in the symmetry project
		if (debug) saveHistory("src/main/java/results/SymmOptimizerHistory.csv");
		return msa;
	}
	
	private void initialize(AFPChain afpChain, Atom[] ca1) throws StructureException, StructureAlignmentException {
		
		//Initialize member variables
		seedAFP = afpChain;
		ca = ca1;
		order = afpChain.getBlockNum();
		subunitLen = afpChain.getOptLen()[0];
		gaps = 0;
		
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
		
		
		//Set the scores and RMSD of the initial state (seed alignment)
		switch (type){
		case CLOSED: 
			updateClosedScore();
			calculatePenaltyDistance();
			updateClosedScore();
			break;
		case OPEN:
			updateOpenScore();
			calculatePenaltyDistance();
			updateOpenScore();
			break;
		}
	}
	
	
	/**
	 *  Optimization method based in a Monte-Carlo approach. Starting from the refined afpChain uses 5 types of moves:
	 *  
	 *  	1- Move Residue: move free residues from one subunit to another.
	 *  	2- Shift Row: if there are enough freePool residues available.
	 *  	3- Expand Block: if there are enough freePool residues available.
	 *  	4- Shrink Block: move a block column to the freePool.
	 *  	5- Split and Shrink Block: split a block in the middle and shrink one column.
	 *  	6- Shift All: shift all the subunits by one position to the right or left.
	 *  	7- Insert gap: insert a gap in a random position and shift the needed residues to right or left.
	 *  
	 * @throws StructureAlignmentException 
	 */
	private void optimizeMC(int maxIter) throws StructureException, StructureAlignmentException{
		
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
			int lastGaps = gaps;
			
			boolean moved = false;
			
			while (!moved){
				//Randomly select one of the steps to modify the alignment
				int move = rnd.nextInt(6);
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
				case 5: moved = insertGap();
						if (debug) System.out.println("did insert gap");
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
					gaps = lastGaps;
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
		generateMultipleAlignment();
	}

	/**
	 * Once the optimization has ended, this method translates the internal data structures to a MultipleAlignment.
	 * @throws StructureAlignmentException 
	 */
	private void generateMultipleAlignment() throws StructureAlignmentException {
		
		//Organize information
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();
		List<Matrix4d> transformations = new ArrayList<Matrix4d>();
		List<String> structureNames = new ArrayList<String>();
		
		for (int i=0; i<order; i++){
			atomArrays.add(ca);
			structureNames.add(seedAFP.getName1());
			Matrix4d transformTimes = new Matrix4d(transformation);
			for (int j=0; j<i; j++){
				transformTimes.mul(transformation);
			}
			transformations.add(transformTimes);
		}
		
		//Initialize the MultipleAlignment to return
		msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setStructureNames(structureNames);
		msa.getEnsemble().setAtomArrays(atomArrays);
		msa.setTransformations(transformations);
		
		//Override the MultipleAlignment with the optimized alignment to return
		msa.setBlockSets(new ArrayList<BlockSet>());
		BlockSet bs = new BlockSetImpl(msa);
		
		//Translate the blocks into a block set
		for (int bkNr=0;bkNr<order;bkNr++){
			List<List<Integer>> newAlgnRes = new ArrayList<List<Integer>>();
			for (int su=0; su<order; su++){
				List<Integer> chain = new ArrayList<Integer>();
				for (int k=0; k<subunitLen; k++){
					chain.add(block.get((su+bkNr)%order).get(k));
				}
				newAlgnRes.add(chain);
			}
			Block bk = new BlockImpl(bs);
			bk.setAlignRes(newAlgnRes);
		}
		
		//Set the scores
		msa.putScore("MC-Score", mcScore);
		msa.putScore("Symm-Score", tmScore);
		msa.putScore("RMSD", rmsd);
		
		//Set the algorithm information
		msa.getEnsemble().setAlgorithmName(seedAFP.getAlgorithmName());
		msa.getEnsemble().setVersion(seedAFP.getVersion());
	}
	
	/**
	 * Method that loops through all the alignment columns and checks that there are no more gaps than the 
	 * maximum allowed.
	 * There must be at least Rmin residues different than null in every alignment column.
	 * In case there is a column with more gaps than allowed it will be shrinked (moved to freePool).
	 * @return true if any columns has been shrinked and false otherwise
	 */
	private boolean checkGaps(){
		
		List<Integer> shrinkColumns = new ArrayList<Integer>();
		//Loop for each column
		for (int res=0; res<subunitLen; res++){
			int gapCount = 0;
			//Loop for each subunit and count the gaps
			for (int su=0; su<order; su++){
				if (block.get(su).get(res) == null) gapCount++;
			}
			if ((order-gapCount)<Rmin){
				//Add the column to the shrink list
				shrinkColumns.add(res);
			}
		}
		
		//Shrink the columns that have more gaps than allowed
		for (int col:shrinkColumns){
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(col);
				block.get(su).remove(col);
				if (residue != null) freePool.get(su).add(residue);
				else gaps--;
				Collections.sort(freePool.get(su));
			}
			subunitLen--;
		}
		
		if (shrinkColumns.size()!=0) return true;
		else return false;
	}
	
	/**
	 * Insert a gap in one of the subunits into a random position in the alignment.
	 * A gap is a null in the block.
	 */
	private boolean insertGap() {

		int su = rnd.nextInt(order); //Select randomly the subunit
		int res = rnd.nextInt(subunitLen); //Position to insert a gap
			
		//Insert the gap at the position
		Integer residueL = block.get(su).get(res);
		if (residueL != null){
			freePool.get(su).add(residueL);
			Collections.sort(freePool.get(su));
			gaps++;
		}
		else return false;  //If there was a gap already in the position.
		
		block.get(su).set(res,null);
		checkGaps();
		return true;
	}

	/**
	 *  Move all the block residues of one subunit one position to the left or right and move the corresponding
	 *  boundary residues from the freePool to the block, and viceversa. 
	 *  The boundaries are determined by any irregularity (either a gap or a discontinuity in the alignment.
	 */
	private boolean shiftRow(){
		
		boolean moved = false;

		int su = rnd.nextInt(order); //Select randomly the subunit that is going to be shifted
		int rl = rnd.nextInt(2);  //Select between moving right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to make the shift
		
		switch(rl){
		case 0: //Move to the right
			//Check if there is a boundary to shidt (a residue different than null to the right)
			while (block.get(su).get(res) == null && res<subunitLen-1) res++;
			if (block.get(su).get(res) == null) return moved;
			
			int leftBoundary = res-1;  //Find the nearest boundary to the left of the pivot
			int leftPrevRes = res;
			while (true){
				if(leftBoundary < 0) break;  //Break if the the left boundary has been found to be the start of the block (=-1)
				else {
					if (block.get(su).get(leftBoundary) == null) break; //Break if there is a gap (this is the boundary)
					else if (block.get(su).get(leftPrevRes) > block.get(su).get(leftBoundary)+1) break;  //Break if there is a discontinuity
				}
				leftPrevRes = leftBoundary;
				leftBoundary--;
			}
			leftBoundary++;
			
			int rightBoundary = res+1;  //Find the nearest boundary to the right of the pivot
			int rightPrevRes = res;
			while (true){
				if(rightBoundary == subunitLen) break;  //Break if the the right boundary has been found (=subunitLen)
				else {
					if (block.get(su).get(rightBoundary) == null) break;  //Break if there is a gap
					else if (block.get(su).get(rightPrevRes)+1 < block.get(su).get(rightBoundary)) break;  //Discontinuity
				}
				rightPrevRes = rightBoundary;
				rightBoundary++;
			}
			rightBoundary--;
			
			//Residues at the boundary
			if (rightBoundary == leftBoundary) return moved;  //If the boundaries are the same they can't be shifted
			Integer residueR0 = block.get(su).get(rightBoundary);
			Integer residueL0 = block.get(su).get(leftBoundary);
			
			//Remove the residue at the right of the block and add it to the freePool
			block.get(su).remove(rightBoundary);
			if (residueR0 != null) freePool.get(su).add(residueR0);
			else throw new IllegalArgumentException("The residue right boundary in shift is null! Cannot be...");
			Collections.sort(freePool.get(su));
			
			//Add the residue at the left of the block from the freePool to the block
			if (residueL0 != null) residueL0 -= 1;
			else throw new IllegalArgumentException("The residue left boundary in shift is null! Cannot be...");
			if (freePool.get(su).contains(residueL0)){
				block.get(su).add(leftBoundary,residueL0);
				freePool.get(su).remove(residueL0);
			}
			else {
				block.get(su).add(leftBoundary,null);
				gaps++;
			}
			
			moved = true;
			break;
			
		case 1: //Move to the left
			//Check if there is a boundary to shidt (a residue different than null to the right)
			while (block.get(su).get(res) == null && res>0) res--;
			if (block.get(su).get(res) == null) return moved;
			
			int leftBoundary1 = res-1;  //Find the nearest boundary to the left of the pivot
			int leftPrevRes1 = res;
			while (true){
				if(leftBoundary1 < 0) break;  //Break if the the left boundary has been found to be the start of the block (=-1)
				else {
					if (block.get(su).get(leftBoundary1) == null) break; //Break if there is a gap (this is the boundary)
					else if (block.get(su).get(leftPrevRes1) > block.get(su).get(leftBoundary1)+1) break;  //Break if there is a discontinuity
				}
				leftPrevRes1 = leftBoundary1;
				leftBoundary1--;
			}
			leftBoundary1++;
			
			int rightBoundary1 = res+1;  //Find the nearest boundary to the right of the pivot
			int rightPrevRes1 = res;
			while (true){
				if(rightBoundary1 == subunitLen) break;  //Break if the the right boundary has been found (=subunitLen)
				else {
					if (block.get(su).get(rightBoundary1) == null) break;  //Break if there is a gap
					else if (block.get(su).get(rightPrevRes1)+1 < block.get(su).get(rightBoundary1)) break;  //Discontinuity
				}
				rightPrevRes1 = rightBoundary1;
				rightBoundary1++;
			}
			rightBoundary1--;
			
			//Residues at the boundary
			if (rightBoundary1 == leftBoundary1) return moved;  //If the boundaries are the same they can't be shifted
			Integer residueR1 = block.get(su).get(rightBoundary1);
			Integer residueL1 = block.get(su).get(leftBoundary1);
			
			//Add the residue at the right of the block from the freePool to the block
			if (residueR1 != null) residueR1 += 1;
			else throw new IllegalArgumentException("The residue right boundary in shift is null! Cannot be...");
			if (freePool.get(su).contains(residueR1)){
				if (rightBoundary1==subunitLen-1) block.get(su).add(residueR1);
				else block.get(su).add(rightBoundary1,residueR1);
				freePool.get(su).remove(residueR1);
			}
			else {
				block.get(su).add(rightBoundary1,null);
				gaps++;
			}
			
			//Remove the residue at the left of the block and add it to the freePool
			block.get(su).remove(leftBoundary1);
			if (residueL1 != null) freePool.get(su).add(residueL1);
			else throw new IllegalArgumentException("The residue left boundary in shift is null! Cannot be...");
			Collections.sort(freePool.get(su));
			
			moved = true;
			break;
		}
		checkGaps();
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
			
			//Loop through all subunits
			for (int su=0; su<order; su++){
				Integer residueEnd = block.get(su).get(subunitLen-1);
				Integer residueStart = block.get(su).get(0);
				//Add a residue to the end
				if (residueEnd == null){
					block.get(su).add(null);
					gaps++;
				} else if (freePool.get(su).contains(residueEnd+1)){
					Integer residueAdd = residueEnd+1;
					block.get(su).add(residueAdd);
					freePool.get(su).remove(residueAdd);
				}
				else {
					block.get(su).add(null);
					gaps++;
				}
				//Delete the residue at the start
				block.get(su).remove(0);
				if (residueStart != null) freePool.get(su).add(residueStart);
				else gaps--;
			}
			moved = true;
			break;
			
		case 1: //Move subunits to the left
			
			//Loop through all subunits
			for (int su=0; su<order; su++){
				Integer residueEnd = block.get(su).get(subunitLen-1);
				Integer residueStart = block.get(su).get(0);
				//Remove the residue to the end
				block.get(su).remove(subunitLen-1);
				if (residueEnd != null) freePool.get(su).add(residueEnd);
				else gaps--;
				
				//Add the residue at the start
				if (residueStart == null){
					block.get(su).add(0,null);
					gaps++;
				} else if (freePool.get(su).contains(residueStart-1)){
					Integer residueAdd = residueStart-1;
					block.get(su).add(0,residueAdd);
					freePool.get(su).remove(residueAdd);
				}
				else {
					block.get(su).add(0,null);
					gaps++;
				}
			}
			moved = true;
			break;
		}
		checkGaps();
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
			int lastIndex = subunitLen-1;
			while (block.get(su).get(lastIndex) == null && lastIndex>0) lastIndex--;
			if (block.get(su).get(lastIndex) == null) return;
			//Check if the residue can be moved, there must exist residues between groups
			if (freePool.get(su).get(n-1) > block.get(su).get(lastIndex)){
				freePool.get(su+1).add(0,freePool.get(su).get(n-1));
				freePool.get(su).remove(n-1);
			}
			break;
			
		case 1:
			//Check if the residue can be moved, there must exist residues between groups
			if (freePool.get(su).size() == 0) return;  //Don't do anything if empty freePool
			int firstIndex = 0;
			while (block.get(su).get(firstIndex) == null && firstIndex<subunitLen-1) firstIndex++;
			if (block.get(su).get(firstIndex) == null) return;
			if (freePool.get(su).get(0) < block.get(su).get(firstIndex)){
				freePool.get(su-1).add(freePool.get(su).get(0));
				freePool.get(su).remove(0);
			}
			break;
		}
	}
	
	/**
	 *  It extends at the beginning or end a group of consecutive aligned residues by moving the residues from the
	 *  freePool to the block. If there are not enough residues in the freePool it introduces gaps.
	 */
	private boolean expandBlock(){
		
		boolean moved = false;
			
		int rl = rnd.nextInt(2);  //Select between expanding right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to expand the subunits
		
		switch (rl) {
		case 0:
			
			//Find the common pivot without any gap to the right of the original pivot
			for (int su=0; su<order; su++){
				while (block.get(su).get(res) == null && res<subunitLen-1) res++;
			}
			
			int rightBoundary = res;
			if (res != subunitLen-1) rightBoundary = findRightBoundary(res);
			
			//Expand the block with the residues at the subunit boundaries
			for (int su=0; su<order; su++){
				Integer residueR = block.get(su).get(rightBoundary);
				if (residueR == null){
					block.get(su).add(null);
					gaps++;
				} else if (freePool.get(su).contains(residueR+1)){
					Integer residueAdd = residueR+1;
					if (rightBoundary == subunitLen-1) block.get(su).add(residueAdd); 
					else block.get(su).add(rightBoundary+1,residueAdd);
					freePool.get(su).remove(residueAdd);
				}
				else {
					block.get(su).add(null);
					gaps++;
				}
			}
			subunitLen++;
			moved = true;
			break;
			
		case 1:
			
			//Find the common pivot without any gap to the left of the original pivot
			for (int su=0; su<order; su++){
				while (block.get(su).get(res) == null && res>0) res--;
			}
			
			int leftBoundary = res;
			if (res!=0) leftBoundary = findLeftBoundary(res);
			
			//Expand the block with the residues at the subunit boundaries
			for (int su=0; su<order; su++){
				Integer residueL = block.get(su).get(leftBoundary);
				if (residueL == null){
					block.get(su).add(leftBoundary,null);
					gaps++;
				} else if (freePool.get(su).contains(residueL+1)){
					Integer residueAdd = residueL-1;
					block.get(su).add(leftBoundary,residueAdd);
					freePool.get(su).remove(residueAdd);
				}
				else {
					block.get(su).add(null);
					gaps++;
				}
			}
			subunitLen++;
			moved = true;
			break;
		}
		if (moved) return !checkGaps();
		return moved;
	}
	
	/**
	 *  It moves aligned residues from the start or end of a group of consecutive aligned residues 
	 *  from the block to the freePool. Shrinks one column of aligned residues at the extremes.
	 */
	private boolean shrinkBlock(){
		
		boolean moved = false;
		if (subunitLen <= Lmin) return moved; //Do not let shrink moves if the subunit is smaller than the minimum length
		
		int rl = rnd.nextInt(2);  //Select between shrink right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to shrink the subunits
		
		switch (rl) {
		case 0:
			
			//Check that there is a boundary
			for (int su=0; su<order; su++){
				while (block.get(su).get(res) == null && res<subunitLen-1) res++;
				if (block.get(su).get(res) == null) return moved;
			}
			int rightBoundary = findRightBoundary(res);
			
			//Shrink the block and add the residues to the freePool
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(rightBoundary);
				block.get(su).remove(rightBoundary);
				if (residue != null) freePool.get(su).add(residue);
				else gaps--;
				Collections.sort(freePool.get(su));
			}
			subunitLen--;
			moved = true;
			break;
			
		case 1:
			
			//Check that there is a boundary
			for (int su=0; su<order; su++){
				while (block.get(su).get(res) == null && res>0) res--;
				if (block.get(su).get(res) == null) return moved;
			}
			int leftBoundary = findLeftBoundary(res);
			
			//Shrink the block and add the residues to the freePool
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(leftBoundary);
				block.get(su).remove(leftBoundary);
				if (residue != null) freePool.get(su).add(residue);
				else gaps--;
				Collections.sort(freePool.get(su));
			}
			subunitLen--;
			moved = true;
			break;
		}
		return moved;
	}
	
	private boolean splitBlock(){
		
		if (subunitLen <= Lmin) return false; //Let split moves everywhere if the subunit is larger than the minimum length
		
		int res = rnd.nextInt(subunitLen); //Residue position to split the subunits
		
		for (int su=0; su<order; su++){
			Integer residue = block.get(su).get(res);
			block.get(su).remove(res);
			if (residue != null) freePool.get(su).add(residue);
			else gaps--;
			Collections.sort(freePool.get(su));
		}
		subunitLen--;
		return true;
	}
	
	/**
	 * Finds the right boundary of consecutive aligned residues in all subunits for a given residue as pivot.
	 */
	private int findRightBoundary(int pivot){
		
		int rightBoundary = pivot+1;  //Find the nearest boundary to the right of the gap
		int rightPrevRes = pivot;
		while (true){
			if(rightBoundary == subunitLen){  //Break if the the right boundary has been found (=subunitLen)
				break;
			} else {
				boolean cont = false;
				boolean brk = true;
				//Loop through all subunits (the discontinuity has to be found in all of them
				for (int su=0; su<order; su++){
					if (block.get(su).get(rightBoundary) == null){
						rightBoundary++;
						cont=true;
						break;
					} else if (block.get(su).get(rightPrevRes)==null){
						brk = false;
						break;
					} else if (block.get(su).get(rightPrevRes)+1 == block.get(su).get(rightBoundary)){
						brk = false;
						break;
					}
				}
				if (cont) continue;
				else if (brk) break;
			}
			rightPrevRes = rightBoundary;
			rightBoundary++;
		}
		return rightBoundary-1;
	}
	
	/**
	 * Fi
	 */
	private int findLeftBoundary(int pivot){
		
		int leftBoundary = pivot-1;  //Find the nearest boundary to the left of the gap
		int leftPrevRes = pivot;
		while (true){
			if(leftBoundary < 0){  //Break if the the left boundary has been found (=0)
				break;
			} else {
				boolean cont = false;
				boolean brk = true;
				for (int su=0; su<order; su++){
					if (block.get(su).get(leftBoundary) == null){
						leftBoundary--;
						cont = true;
						break;
					} else if (block.get(su).get(leftPrevRes)==null){
						brk = false;
						break;
					} else if (block.get(su).get(leftPrevRes) == block.get(su).get(leftBoundary)+1){
						brk = false;
						break;
					}
				}
				if (cont) continue;
				else if (brk) break;
			}
			leftPrevRes = leftBoundary;
			leftBoundary--;
		}
		return leftBoundary+1;
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
		int[] distanceNorm = new int[subunitLen];
		
		//Calculate the aligned atom arrays
		List<Atom> list1 = new ArrayList<Atom>();
		List<Atom> list2 = new ArrayList<Atom>();
		
		for (int j=0; j<order-1; j++){
			for (int k=0; k<subunitLen; k++){
				if (block.get(j).get(k)!=null && block.get((j+1)%order).get(k)!=null){
					list1.add(ca[block.get(j).get(k)]);
					list2.add((Atom) ca[block.get(j+1).get(k)].clone());
				}
			}
		}
		Atom[] arr1 = list1.toArray(new Atom[list1.size()]);
		Atom[] arr2 = list2.toArray(new Atom[list2.size()]);
		
		//Calculate the new transformation information
		SVDSuperimposer svd = new SVDSuperimposer(arr1, arr2);
		transformation = svd.getTransformation();
		
		//Construct all possible transformations of the molecule (order-1) and compare pairwise subunits
		for (int i=0; i<order; i++){
			for (int j=i+1; j<order; j++){
				for (int b=0; b<subunitLen; b++){
					//Consider all the subunit alignment pairs
					if (block.get(i).get(b)!=null && block.get(j).get(b)!=null){
						Atom a1 = ca[block.get(i).get(b)];
						Atom a2 = (Atom) ca[block.get(j).get(b)].clone();
						//Superimpose and generate the distance information
						for (int k=0; k<(j-i); k++) Calc.transform(a2, transformation);
						double distance = Math.abs(Calc.getDistance(a1, a2));
						colDistances[b] += distance;
						distanceNorm[b]++;
						rmsd += distance * distance;
					}
				}
			}
		}
		//Divide the variables for the total number of comparisons to get the average
		int total = 0;
		for (int i=0; i<subunitLen; i++){
			total += distanceNorm[i];
			colDistances[i] /= distanceNorm[i];
		}
		rmsd = Math.sqrt(rmsd/total);
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
		int[] distanceNorm = new int[subunitLen];
		
		//Calculate the aligned atom arrays
		List<Atom> list1 = new ArrayList<Atom>();
		List<Atom> list2 = new ArrayList<Atom>();
		
		for (int j=0; j<order; j++){
			for (int k=0; k<subunitLen; k++){
				if (block.get(j).get(k)!=null && block.get((j+1)%order).get(k)!=null){
					list1.add(ca[block.get(j).get(k)]);
					list2.add((Atom) ca[block.get((j+1)%order).get(k)].clone());
				}
			}
		}
		Atom[] arr1 = list1.toArray(new Atom[list1.size()]);
		Atom[] arr2 = list2.toArray(new Atom[list2.size()]);
		
		//Calculate the new transformation information
		SVDSuperimposer svd = new SVDSuperimposer(arr1, arr2);
		transformation = svd.getTransformation();
		
		//Construct all possible transformations of the molecule (order-1) and compare pairwise subunits
		for (int i=0; i<order; i++){
			for (int j=i+1; j<order; j++){
				for (int b=0; b<subunitLen; b++){
					//Consider all the subunit alignment pairs
					if (block.get(i).get(b)!=null && block.get(j).get(b)!=null){
						Atom a1 = ca[block.get(i).get(b)];
						Atom a2 = (Atom) ca[block.get(j).get(b)].clone();
						//Superimpose and generate the distance information
						for (int k=0; k<(j-i); k++) Calc.transform(a2, transformation);
						double distance = Math.abs(Calc.getDistance(a1, a2));
						colDistances[b] += distance;
						distanceNorm[b]++;
						rmsd += distance * distance;
					}
				}
			}
		}
		//Divide the variables for the total number of comparisons to get the average
		int total = 0;
		for (int i=0; i<subunitLen; i++){
			total += distanceNorm[i];
			colDistances[i] /= distanceNorm[i];
		}
		rmsd = Math.sqrt(rmsd/total);
		mcScore = scoreFunctionCEMC();
		tmScore = scoreFunctionSymmetry();
	}
	
	/**
	 *  Calculates the optimization score from the column average distances.
	 *  Function: sum(M/(d1/d0)^2) minus the penalty A, minus the gap penalty (G) times gaps
	 */
	private double scoreFunctionCEMC(){
		
		double score = 0.0;
		
		//Loop through all the columns
		for (int col=0; col<subunitLen; col++){
			double d1 = colDistances[col];
			double colScore = M/(1+(d1*d1)/(d0*d0));
			colScore -= A;
			colScore *= 2;
			score += colScore;
		}
		return score-gaps*G;
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
						  "4i4q"};  //other
		
		for (String name:names){
			
			System.out.println(name);
			
			AtomCache cache = new AtomCache();
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			CeSymm ceSymm = new CeSymm();
			CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
			params.setRefineMethod(RefineMethod.SINGLE);
			params.setOptimization(false);
			
			AFPChain seedAFP = ceSymm.align(ca1, ca2);
			seedAFP.setName1(name);
			//SingleRefiner refiner = new SingleRefiner();
			//seedAFP = refiner.refine(Arrays.asList(seedAFP), ca1, ca2, 8);
			SymmGapOptimizer optimizer = new SymmGapOptimizer(seedAFP, ca1, SymmetryType.CLOSED, 0);
			MultipleAlignment multAln = optimizer.call();
			
			StructureAlignmentDisplay.display(multAln);
			//SymmetryJmol jmol = new SymmetryJmol(seedAFP, ca1);
			//jmol.setTitle(name);
		}
		
		System.out.println("Finished Alaysis!");
	}
}