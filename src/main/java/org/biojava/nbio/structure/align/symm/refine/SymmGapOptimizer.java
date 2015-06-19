package org.biojava.nbio.structure.align.symm.refine;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Optimizes an alignment by a Monte Carlo score optimization of the subunit multiple alignment.<p>
 * Implements Callable in order to parallelize multiple optimizations.<p>
 * Because gaps are allowed in the subunits, a {@link MultipleAlignment} format is returned.
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
	private static final int Lmin = 15;   //Minimum subunit length
	private int iterFactor = 1000; //Factor to control the max number of iterations of optimization
	private double C = 20; //Probability function constant (probability of acceptance for bad moves)
	
	//Score function parameters
	private static final double M = 20.0; //Maximum score of a match
	private static final double G = 10.0; //Penalty for gaps
	private final double d0 = 5; //Maximum distance that is not penalized
	
	//Alignment Information
	private MultipleAlignment msa;
	private AFPChain seedAFP;
	private Atom[] ca;
	private int order;
	private int subunitLen;
	private int gaps;    //the number of gaps in the alignment
	
	//Multiple Alignment Residues
	private List<List<Integer>> block;     //List to store the residues aligned, in the block. Dimensions are: [order][subunitLen]
	private List<Integer> freePool; 	//List to store the residues not aligned, for rapid check if a residue is aligned.
	
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
	 * Class Constructor. Initializes all the variables needed for the optimization.
	 * @param seedAFP AFPChain with the symmetry subunits split in blocks.
	 * @param ca1
	 * @param type
	 * @param seed
	 * @throws RefinerFailedException 
	 * @throws StructureException 
	 */
	public SymmGapOptimizer(AFPChain seedAFP, Atom[] ca1, SymmetryType type, long seed) throws RefinerFailedException, StructureException {
		
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
		
		optimizeMC(iterFactor*ca.length*order);
		
		//Save the history to the results folder in the symmetry project
		if (debug) saveHistory("src/main/java/results/SymmOptimizerHistory.csv");
		return msa;
	}
	
	private void initialize(AFPChain afpChain, Atom[] ca1) throws StructureException {
		
		//Initialize member variables
		seedAFP = afpChain;
		ca = ca1;
		order = afpChain.getBlockNum();
		subunitLen = afpChain.getOptLen()[0];
		gaps = 0;
		
		//Initialize alignment variables
		block = new ArrayList<List<Integer>>();
		freePool = new ArrayList<Integer>();
		
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
		}
		
		//Add any residue not aligned to the free pool, in the corresponding subunit region (provisional initial state)
		for (int i=0; i<ca.length; i++){
			if (!alreadySeen.contains(i)){
				freePool.add(i);
			}
		}		
		
		//Set the scores and RMSD of the initial state (seed alignment)
		switch (type){
		case CLOSED: 
			updateClosedScore();
			break;
		default:    //case OPEN:
			updateOpenScore();
			break;
		}
	}
	
	
	/**
	 *  Optimization method based in a Monte-Carlo approach. Starting from the refined alignment uses 4 types of moves:
	 *  <p>
	 *  	1- Shift Row: if there are enough freePool residues available.<p>
	 *  	2- Expand Block: add another alignment column if there are residues available.<p>
	 *  	3- Shrink Block: move a block column to the freePool.<p>
	 *  	4- Insert gap: insert a gap in a random position of the alignment.
	 *  
	 */
	private void optimizeMC(int maxIter) throws StructureException{
		
		//Initialize the history variables
		subunitLenHistory = new ArrayList<Integer>();
		rmsdHistory = new ArrayList<Double>();
		scoreHistory = new ArrayList<Double>();
		
		int conv = 0;  //Number of steps without an alignment improvement
		int i = 1;
		int stepsToConverge = maxIter/50;
		
		while (i<maxIter && conv<stepsToConverge){
			
			//Save the state of the system in case the modifications are not favorable
			List<List<Integer>> lastBlock = new ArrayList<List<Integer>>();
			List<Integer> lastFreePool = new ArrayList<Integer>();
			lastFreePool.addAll(freePool);
			for (int k=0; k<order; k++){
				List<Integer> b = new ArrayList<Integer>();
				b.addAll(block.get(k));
				lastBlock.add(b);
			}
			double lastScore = mcScore;
			double lastRMSD = rmsd;
			double lastTMscore = tmScore;
			int lastGaps = gaps;
			
			boolean moved = false;
			
			while (!moved){
				//Randomly select one of the steps to modify the alignment.
				int move = rnd.nextInt(4);
				switch (move){
				case 0: moved = shiftRow();
						if (debug) System.out.println("did shift");
						break;
				case 1: moved = expandBlock();
						if (debug) System.out.println("did expand");
						break;
				case 2:	moved = shrinkBlock();
						if (debug) System.out.println("did shrink");
						break;
				case 3:	moved = insertGap();
						if (debug) System.out.println("did insert gap");
						break;
				}
			}
			
			//Get the properties of the new alignment
			switch (type){
			case CLOSED: 
				updateClosedScore();
				break;
			default:     //case OPEN:
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
	 */
	private void generateMultipleAlignment() {
		
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
		msa.putScore("SymmMC-Score", mcScore);
		msa.putScore("Symm-Score", tmScore);
		msa.putScore(MultipleAlignmentScorer.RMSD, rmsd);
		
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
		
		//Shrink the columns that have more gaps than allowed (from higher indicies to lower ones, to not interfere)
		for (int col=shrinkColumns.size()-1; col>=0; col--){
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(shrinkColumns.get(col));
				block.get(su).remove((int) shrinkColumns.get(col));
				if (residue != null) freePool.add(residue);
				else gaps--;
				Collections.sort(freePool);
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
			freePool.add(residueL);
			Collections.sort(freePool);
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

		int su = rnd.nextInt(order); //Select randomly the subunit that is going to be shifted
		int rl = rnd.nextInt(2);  //Select between moving right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to make the shift
		
		//When the pivot residue is null try to add a residue from the freePool
		if (block.get(su).get(res) == null){
			//Residues not null at the right and left of the pivot null residue
			int rightRes = res;
			int leftRes = res;
			//Find the boundary to the right abd left (a residue different than null)
			while (block.get(su).get(rightRes) == null && rightRes<subunitLen-1) rightRes++;
			while (block.get(su).get(leftRes) == null && leftRes>0) leftRes--;
			
			//If they both are null return because it means that the whole block is null
			if (block.get(su).get(leftRes) == null && block.get(su).get(rightRes) == null) return false;
			else if (block.get(su).get(leftRes) == null){
				//Choose the sequentially previous residue of the known one
				Integer residue = block.get(su).get(rightRes)-1;
				if (freePool.contains(residue)) {
					block.get(su).set(res,residue);
					freePool.remove(residue);
					gaps--;
				} else return false;
			} 
			else if (block.get(su).get(rightRes) == null){
				//Choose the sequentially next residue of the known one
				Integer residue = block.get(su).get(leftRes)+1;
				if (freePool.contains(residue)) {
					block.get(su).set(res,residue);
					freePool.remove(residue);
					gaps--;
				} else return false;
			} 
			else { 
				//If the boundaries are consecutive residues no residue can be added
				if (block.get(su).get(rightRes) == block.get(su).get(leftRes)+1) return false;
				else{
					//Choose randomly a residue in between left and right
					Integer residue = rnd.nextInt(block.get(su).get(rightRes)-block.get(su).get(leftRes)-1) + block.get(su).get(leftRes)+1;
					if (freePool.contains(residue)) {
						block.get(su).set(res,residue);
						freePool.remove(residue);
						gaps--;
					}
					else throw new IllegalStateException("The freePool does not contain a residue in between two aligned positions");
				}
			}
			return true;
		}
		
		//When the residue is different than null try to shift the whole block of consecutive (without gap) residues
		switch(rl){
		case 0: //Move to the right
			
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
			Integer residueR0 = block.get(su).get(rightBoundary);
			Integer residueL0 = block.get(su).get(leftBoundary);
			
			//Remove the residue at the right of the block and add it to the freePool
			block.get(su).remove(rightBoundary);
			if (residueR0 != null) freePool.add(residueR0);
			else throw new IllegalStateException("The residue right boundary in shift is null! Cannot be...");
			Collections.sort(freePool);
			
			//Add the residue at the left of the block from the freePool to the block
			if (residueL0 != null) residueL0 -= 1;
			else throw new IllegalStateException("The residue left boundary in shift is null! Cannot be...");
			if (freePool.contains(residueL0)){
				block.get(su).add(leftBoundary,residueL0);
				freePool.remove(residueL0);
			}
			else {
				block.get(su).add(leftBoundary,null);
				gaps++;
			}
			break;
			
		case 1: //Move to the left
			
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
			Integer residueR1 = block.get(su).get(rightBoundary1);
			Integer residueL1 = block.get(su).get(leftBoundary1);
			
			//Add the residue at the right of the block from the freePool to the block
			if (residueR1 != null) residueR1 += 1;
			else throw new IllegalStateException("The residue right boundary in shift is null! Cannot be...");
			if (freePool.contains(residueR1)){
				if (rightBoundary1==subunitLen-1) block.get(su).add(residueR1);
				else block.get(su).add(rightBoundary1+1,residueR1);
				freePool.remove(residueR1);
			}
			else {
				block.get(su).add(rightBoundary1+1,null);
				gaps++;
			}
			
			//Remove the residue at the left of the block and add it to the freePool
			block.get(su).remove(leftBoundary1);
			if (residueL1 != null) freePool.add(residueL1);
			else throw new IllegalStateException("The residue left boundary in shift is null! Cannot be...");
			Collections.sort(freePool);
			
			break;
		}
		checkGaps();
		return true;
	}
	
	/**
	 *  It extends the alignment one position to the right or to the left of a randomly selected position
	 *  by moving the consecutive residues of each subunit (if present) from the freePool to the block. 
	 *  If there are not enough residues in the freePool it introduces gaps.
	 */
	private boolean expandBlock(){
		
		boolean moved = false;
			
		int rl = rnd.nextInt(2);  //Select between expanding right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot to expand the subunits
		
		switch (rl) {
		case 0:
			
			int rightBoundary = res;
			int[] previousPos = new int[order];
			for (int su=0; su<order; su++) previousPos[su] = -1;
			
			//Search a position to the right that has at minimum Rmin non consecutive residues (otherwise not enough freePool residues to expand)
			while (subunitLen-1>rightBoundary){
				int noncontinuous = 0;
				for (int su=0; su<order; su++){
					if (block.get(su).get(rightBoundary) == null) continue;
					else if (previousPos[su] == -1) previousPos[su] = block.get(su).get(rightBoundary);
					else if (block.get(su).get(rightBoundary) > previousPos[su]+1) noncontinuous++;
				}
				if (noncontinuous < Rmin) rightBoundary++;
				else break;
			}
			rightBoundary--;
			
			//Expand the block with the residues at the subunit boundaries
			for (int su=0; su<order; su++){
				Integer residueR = block.get(su).get(rightBoundary);
				if (residueR == null){
					if (rightBoundary == subunitLen-1) block.get(su).add(null); 
					else block.get(su).add(rightBoundary+1,null);
					gaps++;
				} else if (freePool.contains(residueR+1)){
					Integer residueAdd = residueR+1;
					if (rightBoundary == subunitLen-1) block.get(su).add(residueAdd); 
					else block.get(su).add(rightBoundary+1,residueAdd);
					freePool.remove(residueAdd);
				}
				else {
					if (rightBoundary == subunitLen-1) block.get(su).add(null); 
					else block.get(su).add(rightBoundary+1,null);
					gaps++;
				}
			}
			subunitLen++;
			moved = true;
			break;
			
		case 1:
			
			int leftBoundary = res;
			int[] nextPos = new int[order];
			for (int su=0; su<order; su++) nextPos[su] = -1;
			
			//Search a position to the right that has at minimum Rmin non consecutive residues (otherwise not enough freePool residues to expand)
			while (leftBoundary>0){
				int noncontinuous = 0;
				for (int su=0; su<order; su++){
					if (block.get(su).get(leftBoundary) == null) continue;
					else if (nextPos[su] == -1) nextPos[su] = block.get(su).get(leftBoundary);
					else if (block.get(su).get(leftBoundary) < nextPos[su]-1) noncontinuous++;
				}
				if (noncontinuous < Rmin) leftBoundary--;
				else break;
			}
			
			//Expand the block with the residues at the subunit boundaries
			for (int su=0; su<order; su++){
				Integer residueL = block.get(su).get(leftBoundary);
				if (residueL == null){
					block.get(su).add(leftBoundary,null);
					gaps++;
				} else if (freePool.contains(residueL-1)){
					Integer residueAdd = residueL-1;
					block.get(su).add(leftBoundary,residueAdd);
					freePool.remove(residueAdd);
				}
				else {
					block.get(su).add(leftBoundary, null);
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
	 * Deletes an alignment column at a randomly selected position.
	 */
	private boolean shrinkBlock(){
		
		if (subunitLen <= Lmin) return false; //Let shrink moves only if the subunit is larger than the minimum length
		int res = rnd.nextInt(subunitLen); //Residue position to shrink the subunits
		
		for (int su=0; su<order; su++){
			Integer residue = block.get(su).get(res);
			block.get(su).remove(res);
			if (residue != null) freePool.add(residue);
			else gaps--;
			Collections.sort(freePool);
		}
		subunitLen--;
		checkGaps();
		return true;
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
	 *  Function: sum(M/(d1/d0)^2) minus the gap penalty (G) for each gap
	 */
	private double scoreFunctionCEMC(){
		
		double score = 0.0;
		
		//Loop through all the columns
		for (int col=0; col<subunitLen; col++){
			double d1 = colDistances[col];
			double colScore = M/(1+(d1*d1)/(d0*d0));
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
		return Math.min(Math.max(prob*norm,0.0),1.0);
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
		
		//Easy TIM: "d1i4na_"
		String[] names = { "d1i4na_" };  //Difficult TIMs: "d1hl2a_", "d2fiqa1", "d1eexa_"
		
		for (String name:names){
			
			System.out.println(name);
			
			AtomCache cache = new AtomCache();
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			CeSymm ceSymm = new CeSymm();
			CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
			params.setRefineMethod(RefineMethod.SINGLE);
			params.setSymmetryType(SymmetryType.AUTO);
			params.setOptimization(false);
			
			AFPChain seedAFP = ceSymm.align(ca1, ca2);
			seedAFP.setName1(name);
			//SingleRefiner refiner = new SingleRefiner();
			//seedAFP = refiner.refine(Arrays.asList(seedAFP), ca1, ca2, 8);
			SymmGapOptimizer optimizer = new SymmGapOptimizer(seedAFP, ca1, SymmetryType.CLOSED, 0);
			MultipleAlignment multAln = optimizer.call();
			//SymmOptimizer optimizer = new SymmOptimizer(seedAFP, ca1, SymmetryType.CLOSED, 0);
			//AFPChain multAln = optimizer.call();
			
			MultipleAlignmentDisplay.display(multAln);
			//SymmetryJmol jmol = new SymmetryJmol(multAln, ca1);
			//jmol.setTitle(name);
		}
		
		System.out.println("Finished Alaysis!");
	}
}