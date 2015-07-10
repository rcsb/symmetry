package org.biojava.nbio.structure.align.symm.refine;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentTools;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Optimizes a symmetry alignment by a Monte Carlo score 
 * optimization of the subunit multiple alignment.
 * <p>
 * This algorithm does not use a unfiform distribution for 
 * selecting moves, farther residues have more probability 
 * to be shrinked or gapped.
 * <p>
 * Implements Callable in order to parallelize optimizations.
 * Because gaps are allowed in the subunits, a 
 * {@link MultipleAlignment} format is returned.
 * 
 * @author Aleix Lafita
 * 
 */
public class SymmOptimizer implements Callable<MultipleAlignment> {

	private static final boolean debug = false; //save history and print info
	private SymmetryType type;
	private Random rnd;

	//Optimization parameters
	private static final int Rmin = 2; //min aligned subunits per column
	private static final int Lmin = 15; //min subunit length
	private int maxIterFactor = 100; //max iterations constant
	private double C = 20; //probability of accept bad moves constant

	//Score function parameters
	private static final double Gopen = 20.0; //Penalty for opening gap
	private static final double Gextend = 10.0; //Penalty for extending gaps

	//Alignment Information
	private MultipleAlignment msa;
	private Atom[] atoms;
	private int order;
	private int subunitLen;

	//Multiple Alignment Residues
	private List<List<Integer>> block; //residues aligned
	private List<Integer> freePool; //residues not aligned

	//Optimization information
	private double mcScore; //alignment score to optimize
	private Matrix4d transformation; //elementary symmetry operation

	//Optimization history
	private List<Integer> subunitLenHistory;
	private List<Double> rmsdHistory;
	private List<Double> scoreHistory;

	/**
	 * Constructor with an AFPChain storing a refined symmetry alignment
	 * Initializes all the variables needed for the optimization.
	 * 
	 * @param seedAFP AFPChain with the symmetry subunits split in blocks.
	 * @param atoms
	 * @param type
	 * @param seed
	 * @throws RefinerFailedException 
	 * @throws StructureException 
	 */
	public SymmOptimizer(AFPChain seedAFP, Atom[] atoms, SymmetryType type,
			long seed) throws RefinerFailedException, StructureException {

		//No multiple alignment can be generated if there is only one subunit
		this.order = seedAFP.getBlockNum();
		if (order == 1) {
			throw new RefinerFailedException(
					"Optimization: Non-Symmetric Seed Alignment.");
		}

		this.type = type;
		rnd = new Random(seed);

		initialize(seedAFP, atoms);
	}

	@Override
	public MultipleAlignment call() throws Exception {

		optimizeMC(maxIterFactor*atoms.length);

		//Save the history to the results folder in the symmetry project
		try {
			if (debug) 
				saveHistory("src/main/java/results/SymmOptimizerHistory.csv");
		} catch(FileNotFoundException e) {}

		return msa;
	}

	private void initialize(AFPChain seedAFP, Atom[] ca1) 
			throws StructureException, RefinerFailedException {

		atoms = ca1;
		order = seedAFP.getBlockNum();
		C = 20*order;

		subunitLen = seedAFP.getOptLen()[0];
		if (subunitLen < 1) {
			throw new RefinerFailedException(
					"Optimization: Empty seed alignment!");
		}

		//Initialize MultipleAlignment
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();
		List<String> structureNames = new ArrayList<String>();
		for (int i=0; i<order; i++){
			atomArrays.add(atoms);
			structureNames.add(seedAFP.getName1());
		}
		msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setStructureNames(structureNames);
		msa.getEnsemble().setAtomArrays(atomArrays);
		msa.getEnsemble().setAlgorithmName(seedAFP.getAlgorithmName());
		msa.getEnsemble().setVersion(seedAFP.getVersion());

		//Initialize alignment variables
		block = new ArrayList<List<Integer>>();
		freePool = new ArrayList<Integer>();

		//Store the residues that have been added to the block
		List<Integer> alreadySeen = new ArrayList<Integer>();

		//Generate the initial state of the system
		for (int i=0; i<order; i++){
			List<Integer> residues = new ArrayList<Integer>();
			for (int j=0; j<subunitLen; j++){
				Integer residue = seedAFP.getOptAln()[i][0][j];
				residues.add(residue);
				alreadySeen.add(residue);
			}
			block.add(residues);
		}

		//Add any residue not aligned to the free pool
		for (int i=0; i<atoms.length; i++){
			if (!alreadySeen.contains(i)){
				freePool.add(i);
			}
		}		

		//Set the MC score and RMSD of the initial state (seed alignment)
		updateTransformation();
		subunitMultipleAlignment();
		mcScore = MultipleAlignmentScorer.getMCScore(msa, Gopen, Gextend);
	}

	/**
	 *  Optimization method based in a Monte-Carlo approach. 
	 *  Starting from the refined alignment uses 4 types of moves:
	 *  <p>
	 *  	1- Shift Row: if there are enough freePool residues available.<p>
	 *  	2- Expand Block: add another alignment column if there are residues
	 *  		available.<p>
	 *  	3- Shrink Block: move a block column to the freePool.<p>
	 *  	4- Insert gap: insert a gap in a position of the alignment.
	 *  
	 */
	private void optimizeMC(int maxIter) throws StructureException{

		//Initialize the history variables
		subunitLenHistory = new ArrayList<Integer>();
		rmsdHistory = new ArrayList<Double>();
		scoreHistory = new ArrayList<Double>();

		int conv = 0;  //Number of steps without an alignment improvement
		int i = 1;
		int stepsToConverge = Math.max(maxIter/50,1000);

		while (i<maxIter && conv<stepsToConverge){

			//Save the state of the system
			List<List<Integer>> lastBlock = new ArrayList<List<Integer>>();
			List<Integer> lastFreePool = new ArrayList<Integer>();
			lastFreePool.addAll(freePool);
			for (int k=0; k<order; k++){
				List<Integer> b = new ArrayList<Integer>();
				b.addAll(block.get(k));
				lastBlock.add(b);
			}
			double lastScore = mcScore;

			boolean moved = false;

			while (!moved){
				//Randomly select one of the steps to modify the alignment. 
				//Because of biased moves, the probabilities are not the same
				double move = rnd.nextDouble();
				if (move < 0.4){
					moved = shiftRow();
					if (debug) System.out.println("did shift");
				}
				else if (move < 0.7){
					moved = expandBlock();
					if (debug) System.out.println("did expand");
				}
				else if (move < 0.85){
					moved = shrinkBlock();
					if (debug) System.out.println("did shrink");
				}
				else {
					moved = insertGap();
					if (debug) System.out.println("did insert gap");
				}
			}

			//Get the properties of the new alignment
			updateTransformation();
			subunitMultipleAlignment();
			mcScore = MultipleAlignmentScorer.getMCScore(msa, Gopen, Gextend);

			//Calculate change in the optimization Score
			double AS = mcScore-lastScore;
			double prob=1.0;

			if (AS<0){

				//Probability of accepting bad move
				prob = probabilityFunction(AS,i,maxIter);
				double p = rnd.nextDouble();

				//Reject the move
				if (p>prob){
					block = lastBlock;
					freePool = lastFreePool;
					subunitLen = block.get(0).size();
					mcScore = lastScore;
					conv ++; //no change in score if rejected

				} else conv = 0; //if accepted

			} else conv=0; //if positive change

			if (debug) {
				System.out.println(i+": --prob: "+prob+", --score: "+AS+
						", --conv: "+conv);

				if (i%100==1){
					//Get the correct superposition again
					updateTransformation();
					subunitMultipleAlignment();
					double rmsd = MultipleAlignmentScorer.getRMSD(msa);

					subunitLenHistory.add(subunitLen);
					rmsdHistory.add(rmsd);
					scoreHistory.add(mcScore);
				}	
			}

			i++;
		}
		//Superimpose and calculate scores
		updateTransformation();
		subunitMultipleAlignment();
		mcScore = MultipleAlignmentScorer.getMCScore(msa, Gopen, Gextend);
		double tmScore = MultipleAlignmentScorer.getAvgTMScore(msa) * order;
		double rmsd = MultipleAlignmentScorer.getRMSD(msa);

		//Set the scores
		msa.putScore(MultipleAlignmentScorer.MC_SCORE, mcScore);
		msa.putScore(MultipleAlignmentScorer.AVGTM_SCORE, tmScore);
		msa.putScore(MultipleAlignmentScorer.RMSD, rmsd);
	}

	/**
	 * This method translates the internal data structures to a 
	 * MultipleAlignment of the subunits only in order to use 
	 * the normal methods to score MultipleAlignments.
	 */
	private void subunitMultipleAlignment() {

		msa.clear();
		//Create transformations from the symmetry operation
		List<Matrix4d> transformations = new ArrayList<Matrix4d>();		
		for (int i=0; i<order; i++){
			Matrix4d transformTimes = new Matrix4d();
			transformTimes.setIdentity();
			for (int j=0; j<i; j++) transformTimes.mul(transformation);
			transformations.add(transformTimes);
		}
		msa.setTransformations(transformations);

		//Override the alignment with the optimized alignment of the subunits
		msa.setBlockSets(new ArrayList<BlockSet>());
		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(block);
	}

	/**
	 * Method that loops through all the alignment columns and checks 
	 * that there are no more gaps than the maximum allowed: Rmin.
	 * <p>
	 * There must be at least Rmin residues different than null in 
	 * every alignment column.In case there is a column with more 
	 * gaps than allowed it will be shrinked (moved to freePool).
	 * 
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
		for (int col=shrinkColumns.size()-1; col>=0; col--){
			for (int su=0; su<order; su++){
				Integer residue = block.get(su).get(shrinkColumns.get(col));
				block.get(su).remove((int) shrinkColumns.get(col));
				if (residue != null) freePool.add(residue);
				Collections.sort(freePool);
			}
			subunitLen--;
		}

		if (shrinkColumns.size()!=0) return true;
		else return false;
	}

	/**
	 * Insert a gap in one of the subunits into selected position 
	 * (by higher distances) in the alignment. Calculates the average
	 * residue distance to make the choice.
	 * A gap is a null in the block.
	 * 
	 * @throws StructureException 
	 */
	private boolean insertGap() throws StructureException {

		//Let gaps only if the subunit is larger than the minimum length
		if (subunitLen <= Lmin) return false;

		//Select residue by maximum distance
		updateTransformation();
		subunitMultipleAlignment();
		Matrix residueDistances = 
				MultipleAlignmentTools.getAverageResidueDistances(msa);

		double maxDist = Double.MIN_VALUE;
		int su = 0;
		int res = 0;
		for (int col=0; col<subunitLen; col++){
			for (int s=0; s<order; s++){
				if (residueDistances.get(s, col) != -1){
					if (residueDistances.get(s, col) > maxDist){
						//geometric distribution
						if (rnd.nextDouble() > 0.5) {
							su = s;
							res = col;
							maxDist = residueDistances.get(s, col);
						}
					}
				}
			}
		}

		//Insert the gap at the position
		Integer residueL = block.get(su).get(res);
		if (residueL != null){
			freePool.add(residueL);
			Collections.sort(freePool);
		}
		else return false;  //If there was a gap already in the position.

		block.get(su).set(res,null);
		checkGaps();
		return true;
	}

	/**
	 *  Move all the block residues of one subunit one position 
	 *  to the left or right and move the corresponding boundary 
	 *  residues from the freePool to the block, and viceversa.
	 *  <p>
	 *  The boundaries are determined by any irregularity 
	 *  (either a gap or a discontinuity in the alignment.
	 */
	private boolean shiftRow(){

		int su = rnd.nextInt(order); //Select the subunit
		int rl = rnd.nextInt(2);  //Select between moving right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot

		//When the pivot residue is null try to add a residue from the freePool
		if (block.get(su).get(res) == null){

			int right = res;
			int left = res;
			//Find the boundary to the right abd left
			while (block.get(su).get(right) == null && right<subunitLen-1) {
				right++;
			}
			while (block.get(su).get(left) == null && left>0) {
				left--;
			}

			//If they both are null the whole block is null
			if (block.get(su).get(left) == null && 
					block.get(su).get(right) == null) {
				return false;
			}
			else if (block.get(su).get(left) == null){
				//Choose the sequentially previous residue of the known one
				Integer residue = block.get(su).get(right)-1;
				if (freePool.contains(residue)) {
					block.get(su).set(res,residue);
					freePool.remove(residue);
				} else return false;
			} 
			else if (block.get(su).get(right) == null){
				//Choose the sequentially next residue of the known one
				Integer residue = block.get(su).get(left)+1;
				if (freePool.contains(residue)) {
					block.get(su).set(res,residue);
					freePool.remove(residue);
				} else return false;
			} 
			else { 
				//If boundaries are consecutive swap null and position (R or L)
				if (block.get(su).get(right) == block.get(su).get(left)+1){
					switch(rl){
					case 0: //to the right
						block.get(su).set(right-1,block.get(su).get(right));
						block.get(su).set(right, null);
						break;
					case 1: //to the left
						block.get(su).set(left+1,block.get(su).get(left));
						block.get(su).set(left, null);
						break;
					}
				}
				else{
					//Choose randomly a residue in between left and right to add
					Integer residue = rnd.nextInt(block.get(su).get(right)-
									block.get(su).get(left)-1) + 
									block.get(su).get(left)+1;
					
					if (freePool.contains(residue)) {
						block.get(su).set(res,residue);
						freePool.remove(residue);
					}
				}
			}
			return true;
		}

		//When the residue is different than null
		switch(rl){
		case 0: //Move to the right

			int leftBoundary = res-1;
			int leftPrevRes = res;
			while (true){
				if(leftBoundary < 0) break;
				else {
					if (block.get(su).get(leftBoundary) == null) {
						break; //gap
					}
					else if (block.get(su).get(leftPrevRes) 
							> block.get(su).get(leftBoundary)+1) {
						break; //discontinuity
					}
				}
				leftPrevRes = leftBoundary;
				leftBoundary--;
			}
			leftBoundary++;

			int rightBoundary = res+1;
			int rightPrevRes = res;
			while (true){
				if(rightBoundary == subunitLen) break;
				else {
					if (block.get(su).get(rightBoundary) == null) {
						break;  //gap
					}
					else if (block.get(su).get(rightPrevRes)+1 
							< block.get(su).get(rightBoundary)){ 
						break;  //discontinuity
					}
				}
				rightPrevRes = rightBoundary;
				rightBoundary++;
			}
			rightBoundary--;

			//Residues at the boundary
			Integer residueR0 = block.get(su).get(rightBoundary);
			Integer residueL0 = block.get(su).get(leftBoundary);

			//Remove residue at the right of the block and add to the freePool
			block.get(su).remove(rightBoundary);
			if (residueR0 != null) {
				freePool.add(residueR0);
				Collections.sort(freePool);
			}

			//Add the residue at the left of the block
			residueL0 -= 1; //cannot be null, throw exception if it is
			if (freePool.contains(residueL0)){
				block.get(su).add(leftBoundary,residueL0);
				freePool.remove(residueL0);
			}
			else {
				block.get(su).add(leftBoundary,null);
			}
			break;

		case 1: //Move to the left

			int leftBoundary1 = res-1;
			int leftPrevRes1 = res;
			while (true){
				if(leftBoundary1 < 0) break;
				else {
					if (block.get(su).get(leftBoundary1) == null) {
						break; //gap
					}
					else if (block.get(su).get(leftPrevRes1) 
							> block.get(su).get(leftBoundary1)+1) {
						break;  //discontinuity
					}
				}
				leftPrevRes1 = leftBoundary1;
				leftBoundary1--;
			}
			leftBoundary1++;

			int rightBoundary1 = res+1;
			int rightPrevRes1 = res;
			while (true){
				if(rightBoundary1 == subunitLen) break;
				else {
					if (block.get(su).get(rightBoundary1) == null) {
						break;  //gap
					}
					else if (block.get(su).get(rightPrevRes1)+1 
							< block.get(su).get(rightBoundary1)) {
						break;  //discontinuity
					}
				}
				rightPrevRes1 = rightBoundary1;
				rightBoundary1++;
			}
			rightBoundary1--;

			//Residues at the boundary
			Integer residueR1 = block.get(su).get(rightBoundary1);
			Integer residueL1 = block.get(su).get(leftBoundary1);

			//Add the residue at the right of the block
			residueR1 += 1; //cannot be null
			if (freePool.contains(residueR1)){
				if (rightBoundary1==subunitLen-1) block.get(su).add(residueR1);
				else block.get(su).add(rightBoundary1+1,residueR1);
				freePool.remove(residueR1);
			}
			else {
				block.get(su).add(rightBoundary1+1,null);
			}

			//Remove the residue at the left of the block
			block.get(su).remove(leftBoundary1);
			freePool.add(residueL1);
			Collections.sort(freePool);
			break;
		}
		checkGaps();
		return true;
	}

	/**
	 *  It extends the alignment one position to the right 
	 *  or to the left of a randomly selected position
	 *  by moving the consecutive residues of each subunit 
	 *  (if present) from the freePool to the block.
	 *  <p>
	 *  If there are not enough residues in the freePool 
	 *  it introduces gaps.
	 */
	private boolean expandBlock(){

		boolean moved = false;

		int rl = rnd.nextInt(2);  //Select between right (0) or left (1)
		int res = rnd.nextInt(subunitLen); //Residue as a pivot

		switch (rl) {
		case 0:

			int rightBoundary = res;
			int[] previousPos = new int[order];
			for (int su=0; su<order; su++) previousPos[su] = -1;

			//Search a position to the right that has at minimum Rmin
			while (subunitLen-1>rightBoundary){
				int noncontinuous = 0;
				for (int su=0; su<order; su++){
					if (block.get(su).get(rightBoundary) == null) {
						continue;
					}
					else if (previousPos[su] == -1) {
						previousPos[su] = block.get(su).get(rightBoundary);
					}
					else if (block.get(su).get(rightBoundary) 
							> previousPos[su]+1) {
						noncontinuous++;
					}
				}
				if (noncontinuous < Rmin) rightBoundary++;
				else break;
			}
			if (rightBoundary > 0) rightBoundary--;

			//Expand the block with the residues at the subunit boundaries
			for (int su=0; su<order; su++){
				Integer residueR = block.get(su).get(rightBoundary);
				if (residueR == null){
					if (rightBoundary == subunitLen-1) block.get(su).add(null); 
					else block.get(su).add(rightBoundary+1,null);
				} else if (freePool.contains(residueR+1)){
					Integer residueAdd = residueR+1;
					if (rightBoundary == subunitLen-1) {
						block.get(su).add(residueAdd); 
					}
					else block.get(su).add(rightBoundary+1,residueAdd);
					freePool.remove(residueAdd);
				}
				else {
					if (rightBoundary == subunitLen-1) block.get(su).add(null);
					else block.get(su).add(rightBoundary+1,null);
				}
			}
			subunitLen++;
			moved = true;
			break;

		case 1:

			int leftBoundary = res;
			int[] nextPos = new int[order];
			for (int su=0; su<order; su++) nextPos[su] = -1;

			//Search a position to the right that has at minimum Rmin
			while (leftBoundary>0){
				int noncontinuous = 0;
				for (int su=0; su<order; su++){
					if (block.get(su).get(leftBoundary) == null){
						continue;
					}
					else if (nextPos[su] == -1){
						nextPos[su] = block.get(su).get(leftBoundary);
					}
					else if (block.get(su).get(leftBoundary) 
							< nextPos[su]-1) {
						noncontinuous++;
					}
				}
				if (noncontinuous < Rmin) leftBoundary--;
				else break;
			}

			//Expand the block with the residues at the subunit boundaries
			for (int su=0; su<order; su++){
				Integer residueL = block.get(su).get(leftBoundary);
				if (residueL == null){
					block.get(su).add(leftBoundary,null);
				} else if (freePool.contains(residueL-1)){
					Integer residueAdd = residueL-1;
					block.get(su).add(leftBoundary,residueAdd);
					freePool.remove(residueAdd);
				}
				else {
					block.get(su).add(leftBoundary, null);
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
	 * @throws StructureException 
	 */
	private boolean shrinkBlock() throws StructureException{

		//Let shrink moves only if the subunit is larger enough
		if (subunitLen <= Lmin) return false;

		//Select column by maximum distance
		updateTransformation();
		subunitMultipleAlignment();
		Matrix residueDistances = 
				MultipleAlignmentTools.getAverageResidueDistances(msa);
		
		double maxDist = Double.MIN_VALUE;
		double[] colDistances = new double[subunitLen];
		int res = 0;
		for (int col=0; col<subunitLen; col++){
			int normalize = 0;
			for (int s=0; s<order; s++){
				if (residueDistances.get(s, col) != -1){
					colDistances[col] += residueDistances.get(s, col);
					normalize++;
				}
			}
			colDistances[col] /= normalize;
			if (colDistances[col] > maxDist){
				//geometric distribution
				if (rnd.nextDouble() > 0.5) {
					maxDist = colDistances[col];
					res = col;
				}
			}
		}

		for (int su=0; su<order; su++){
			Integer residue = block.get(su).get(res);
			block.get(su).remove(res);
			if (residue != null) freePool.add(residue);
			Collections.sort(freePool);
		}
		subunitLen--;
		checkGaps();
		return true;
	}

	/**
	 *  Calculates the symmetry operation Matrix (transformation) 
	 *  for the new alignment.
	 */
	private void updateTransformation() throws StructureException {

		//Calculate the aligned atom arrays
		List<Atom> list1 = new ArrayList<Atom>();
		List<Atom> list2 = new ArrayList<Atom>();

		switch (type) {
		case CLOSED:
			for (int j=0; j<order; j++){
				for (int k=0; k<subunitLen; k++){
					if (block.get(j).get(k)!=null && 
							block.get((j+1)%order).get(k)!=null){
						
						list1.add(atoms[block.get(j).get(k)]);
						list2.add((Atom) 
								atoms[block.get((j+1)%order).get(k)].clone());
					}
				}
			}
			Atom[] arr1 = list1.toArray(new Atom[list1.size()]);
			Atom[] arr2 = list2.toArray(new Atom[list2.size()]);

			//Calculate the new transformation information
			SVDSuperimposer svd = new SVDSuperimposer(arr1, arr2);
			transformation = svd.getTransformation();
			break;

		default: //case OPEN:
			for (int j=0; j<order-1; j++){
				for (int k=0; k<subunitLen; k++){
					if (block.get(j).get(k)!=null && 
							block.get((j+1)%order).get(k)!=null){
						
						list1.add(atoms[block.get(j).get(k)]);
						list2.add((Atom) atoms[block.get(j+1).get(k)].clone());
					}
				}
			}
			Atom[] arr3 = list1.toArray(new Atom[list1.size()]);
			Atom[] arr4 = list2.toArray(new Atom[list2.size()]);

			//Calculate the new transformation information
			SVDSuperimposer svd2 = new SVDSuperimposer(arr3, arr4);
			transformation = svd2.getTransformation();
			break;
		}
	}	

	/**
	 *  Calculates the probability of accepting a bad move 
	 *  given the iteration step and the score change.
	 *  <p>
	 *  Function: p=(C-AS)/(C*sqrt(step))
	 *  Added a normalization factor so that the probability 
	 *  approaches 0 when the maxIter is reached.
	 */
	private double probabilityFunction(double AS, int m, int maxIter) {

		double prob = (C+AS)/(C*Math.sqrt(m));
		double norm = (1-(m*1.0)/maxIter);  //Normalization factor
		return Math.min(Math.max(prob*norm,0.0),1.0);
	}

	/**
	 * Save the evolution of the optimization process as a csv file.
	 */
	private void saveHistory(String filePath) throws IOException {

		FileWriter writer = new FileWriter(filePath);
		writer.append("Step,Length,RMSD,Score\n");

		for (int i=0; i<subunitLenHistory.size(); i++){
			writer.append(
					i*100+","+subunitLenHistory.get(i)+
					","+rmsdHistory.get(i)+","+scoreHistory.get(i)+"\n");
		}

		writer.flush();
		writer.close();
	}

	public static void main(String[] args) throws Exception{

		//Easy TIM: "d1i4na_"
		//Crystallins: 4GCR, d4gcra1, d4gcra2
		//Aspartic proteinases: 
		//Difficult TIMs: "d1hl2a_", "d2fiqa1", "d1eexa_"
		String[] names = {"4i4q"};
		AtomCache cache = new AtomCache();

		for (String name:names){

			System.out.println(name);
			List<Atom[]> atoms = new ArrayList<Atom[]>();
			atoms.add(cache.getAtoms(name));

			CeSymm ceSymm = new CeSymm();

			CESymmParameters params = 
					(CESymmParameters) ceSymm.getParameters();
			
			params.setRefineMethod(RefineMethod.SINGLE);
			params.setSymmetryType(SymmetryType.AUTO);
			params.setOptimization(true);
			//params.setSeed(10);

			MultipleAlignment symmetry = ceSymm.align(atoms);

			new SymmetryJmol(symmetry);
		}
	}
}
