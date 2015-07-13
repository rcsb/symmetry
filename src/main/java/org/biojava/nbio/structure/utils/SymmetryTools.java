package org.biojava.nbio.structure.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.helper.AlignTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Utility methods for the internal symmetry identification and manipulation.
 * <p>
 * Methods include: blank out regions of DP Matrix, build symmetry graphs,
 * get rotation symmetry angles, split subunits in chains, convert between
 * storing formats.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class SymmetryTools {

	//there won't be an instance of this
	private SymmetryTools(){}

	/**
	 * Returns the "reset value" for graying out the main diagonal. 
	 * If we're blanking out the main diagonal, this value is always 
	 * Integer.MIN_VALUE.
	 * <p>
	 * This is possible if {@code gradientPolyCoeff = {Integer.MIN_VALUE}} 
	 * and {@code gradientExpCoeff = 0}.
	 * 
	 * @param unpenalizedScore
	 * @param nResFromMainDiag
	 * @param gradientPolyCoeff
	 * @param gradientExpCoeff
	 * @return
	 */
	private static double getResetVal(double unpenalizedScore, 
			double nResFromMainDiag, double[] gradientPolyCoeff, 
			double gradientExpCoeff) {

		if (unpenalizedScore == Double.NaN) 
			return 0; // what else?

		//We can actually return a positive value if this is high enough
		double updateVal = unpenalizedScore;
		updateVal -= gradientExpCoeff * Math.pow(Math.E, -nResFromMainDiag);
		for (int p = 0; p < gradientPolyCoeff.length; p++) {
			updateVal -= gradientPolyCoeff[gradientPolyCoeff.length-1-p]
					* Math.pow(nResFromMainDiag, -p);
		}
		return updateVal;
	}

	/**
	 * Grays out the main diagonal of a duplicated distance matrix.
	 * 
	 * @param ca2
	 * @param rows Number of rows
	 * @param cols Number of original columns
	 * @param calculator Used to get the matrix if origM is null
	 * @param origM starting matrix. If null, uses 
	 * 			{@link CECalculator#getMatMatrix()}
	 * @param blankWindowSize Width of section to gray out
	 * @param gradientPolyCoeff
	 * @param gradientExpCoeff
	 * @return
	 */
	public static Matrix grayOutCEOrig(Atom[] ca2, int rows, int cols,
			CECalculator calculator, Matrix origM, int blankWindowSize, 
			double[] gradientPolyCoeff, double gradientExpCoeff) {

		if (origM == null){
			origM = new Matrix(calculator.getMatMatrix());
		}

		//symmetry hack, disable main diagonal

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				double resetVal = getResetVal(origM.get(i, j), 
						diff, gradientPolyCoeff, gradientExpCoeff);

				if ( diff < blankWindowSize ){
					origM.set(i,j, origM.get(i, j) + resetVal);

				}
				int diff2 = Math.abs(i-(j-ca2.length/2)); // other side

				double resetVal2 = getResetVal(origM.get(i, j), 
						diff2, gradientPolyCoeff, gradientExpCoeff);

				if ( diff2 < blankWindowSize ){
					origM.set(i,j, origM.get(i, j) + resetVal2);

				}
			}
		}
		return origM;
	}

	public static Matrix grayOutPreviousAlignment(AFPChain afpChain, 
			Atom[] ca2, int rows, int cols, CECalculator calculator, 
			Matrix max, int blankWindowSize, double[] gradientPolyCoeff, 
			double gradientExpCoeff) {

		max =  grayOutCEOrig(ca2, rows, cols, calculator, max,  
				blankWindowSize, gradientPolyCoeff, gradientExpCoeff);

		double[][] dist1 = calculator.getDist1();
		double[][] dist2 = calculator.getDist2();

		int[][][] optAln = afpChain.getOptAln();
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();

		// ca2 is circularly permutated
		int breakPoint = ca2.length / 2;
		for (int bk = 0; bk < blockNum; bk++)       {

			for ( int i=0;i< optLen[bk];i++){
				int pos1 = optAln[bk][0][i];
				int pos2 = optAln[bk][1][i];

				int dist = blankWindowSize/2 ;
				int start1 = Math.max(pos1-dist,0);
				int start2 = Math.max(pos2-dist,0);
				int end1 = Math.min(pos1+dist, rows-1);
				int end2 = Math.min(pos2+dist, cols-1);

				for ( int i1 = start1; i1< end1 ; i1++){

					// blank diagonal of dist1
					for ( int k=0; k < blankWindowSize/2 ; k ++){
						if ( i1-k >= 0) {
							double resetVal = getResetVal(max.get(i1-k, i1-k),
									0, gradientPolyCoeff, gradientExpCoeff);
							dist1[i1-k][i1-k] = resetVal;
						} else if ( i1+k < rows) {
							double resetVal = getResetVal(max.get(i1+k, i1+k),
									0, gradientPolyCoeff, gradientExpCoeff);
							dist1[i1+k][i1+k] = resetVal;
						}

					}

					for ( int j2 = start2 ; j2 < end2 ; j2++){
						double resetVal = getResetVal(max.get(i1, j2), 
								Math.abs(i1-j2), gradientPolyCoeff, 
								gradientExpCoeff);
						max.set(i1,j2,resetVal);
						if ( j2 < breakPoint) {
							double resetVal2 = getResetVal(max.get(i1, 
									j2+breakPoint), 
									Math.abs(i1-(j2+breakPoint)), 
									gradientPolyCoeff, 
									gradientExpCoeff);
							max.set(i1,j2+breakPoint,resetVal2);
						} else {
							double resetVal2 = getResetVal(max.get(i1, 
									j2-breakPoint), 
									Math.abs(i1-(j2-breakPoint)), 
									gradientPolyCoeff, 
									gradientExpCoeff);
							max.set(i1,j2-breakPoint,resetVal2);
						}
						for ( int k=0; k <blankWindowSize/2 ; k ++){
							if ( j2-k >=0) {
								if( j2-k < breakPoint ) {
									double resetVal2 = getResetVal(
											max.get(j2-k, j2-k), 0, 
											gradientPolyCoeff, 
											gradientExpCoeff);
									dist2[j2-k][j2-k] = resetVal2;
								} else {
									double resetVal2 = getResetVal(
											max.get(j2-k-breakPoint, j2-k), 0,
											gradientPolyCoeff, 
											gradientExpCoeff);
									dist2[j2-k-breakPoint][j2-k-breakPoint] = 
											resetVal2;
								}
							} else if ( j2+k < cols) {
								if( j2+k < breakPoint) {
									double resetVal2 = getResetVal(
											max.get(j2+k, j2+k), 0, 
											gradientPolyCoeff, 
											gradientExpCoeff);
									dist2[j2+k][j2+k] = resetVal2;
								} else {
									double resetVal2 = getResetVal(
											max.get(j2+k-breakPoint, j2+k), 0,
											gradientPolyCoeff, 
											gradientExpCoeff);
									dist2[j2+k-breakPoint][j2+k-breakPoint] =
											resetVal2;
								}
							}
						}
					}
				}

			}
		}
		calculator.setDist1(dist1);
		calculator.setDist2(dist2);
		return max;

	}

	public Matrix getDkMatrix(Atom[] ca1, Atom[] ca2,int fragmentLength,
			double[] dist1, double[] dist2, int rows, int cols) {
		
		Matrix diffDistMax =  Matrix.identity(ca1.length, ca2.length);

		for ( int i = 0 ; i< rows; i++){
			double score1 = 0;
			for ( int x=0 ; x < fragmentLength ; x++){
				score1 += dist1[i+x];
			}
			for ( int j = 0 ; j < cols ; j++){
				double score2 = 0;
				for ( int y=0 ; y < fragmentLength ; y++){
					score2 += dist2[j+y];
				}

				// if the intramolecular distances are very similar
				// the two scores should be similar, 
				//i.e. the difference is close to 0
				diffDistMax.set(i,j, Math.abs(score1-score2));
			}
		}


		// symmetry hack, disable main diagonal

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				if ( diff < 15 ){
					diffDistMax.set(i,j, 99);
				}
				int diff2 = Math.abs(i-(j-ca2.length/2));
				if ( diff2 < 15 ){
					diffDistMax.set(i,j, 99);
				}
			}
		}
		return diffDistMax;

	}

	public static Matrix blankOutPreviousAlignment(AFPChain afpChain, 
			Atom[] ca2,	int rows, int cols, CECalculator calculator, 
			Matrix max, int blankWindowSize) {
		return grayOutPreviousAlignment(afpChain, ca2, rows, cols, calculator,
				max, blankWindowSize, new double[] {Integer.MIN_VALUE}, 0.0);
	}

	public static Matrix blankOutCEOrig(Atom[] ca2, int rows, int cols,
			CECalculator calculator, Matrix origM, int blankWindowSize) {
		return grayOutCEOrig(ca2, rows, cols, calculator, origM, 
				blankWindowSize, new double[] {Integer.MIN_VALUE}, 0.0);
	}

	@Deprecated
	public static Atom[] cloneAtoms(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		/*Atom[] ca2clone = new Atom[ca2.length];

		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom

			ca2clone[pos] = g.getAtom(a.getName());

			pos++;
		}

		return ca2clone;*/
		return StructureTools.cloneAtomArray(ca2);
	}

	public static Matrix getDkMatrix(Atom[] ca1, Atom[] ca2, int k, 
			int fragmentLength) {
		
		double[] dist1 = AlignTools.getDiagonalAtK(ca1, k);
		double[] dist2 = AlignTools.getDiagonalAtK(ca2, k);

		int rows = ca1.length - fragmentLength - k + 1;
		int cols = ca2.length - fragmentLength - k + 1;

		// Matrix that tracks similarity of a fragment of length fragmentLength
		// starting a position i,j.

		Matrix m2 = new Matrix(rows,cols); 

		for ( int i = 0 ; i< rows; i++){
			double score1 = 0;
			for ( int x=0 ; x < fragmentLength ; x++){
				score1 += dist1[i+x];
			}
			for ( int j = 0 ; j < cols ; j++){
				double score2 = 0;
				for ( int y=0 ; y < fragmentLength ; y++){
					score2 += dist2[j+y];
				}	

				// if the intramolecular distances are very similar
				// the two scores should be similar, 
				//i.e. the difference is close to 0
				m2.set(i,j, Math.abs(score1-score2));
			}
		}
		return m2;
	}


	public static boolean[][] blankOutBreakFlag(AFPChain afpChain,
			Atom[] ca2, int rows, int cols, CECalculator calculator,
			boolean[][] breakFlag, int blankWindowSize) {


		int[][][] optAln = afpChain.getOptAln();
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();

		// ca2 is circularly permutated at this point.
		int breakPoint = ca2.length / 2;

		for(int bk = 0; bk < blockNum; bk ++)       {

			//Matrix m= afpChain.getBlockRotationMatrix()[bk];
			//Atom shift = afpChain.getBlockShiftVector()[bk];
			for ( int i=0;i< optLen[bk];i++){
				int pos1 = optAln[bk][0][i];
				int pos2 = optAln[bk][1][i];
				// blank out area around these positions...

				int dist = blankWindowSize ;
				int start1 = Math.max(pos1-dist,0);
				int start2 = Math.max(pos2-dist,0);
				int end1 = Math.min(pos1+dist, rows-1);
				int end2 = Math.min(pos2+dist, cols-1);

				for ( int i1 = start1; i1< end1 ; i1++){

					for ( int j2 = start2 ; j2 < end2 ; j2++){
						//System.out.println(i1 + " " + j2 + " (***)");
						breakFlag[i1][j2] = true;
						if ( j2 < breakPoint) {
							breakFlag[i1][j2+ breakPoint ] = true;
						}
					}
				}

			}
		}

		return breakFlag;
	}

	/**
	 * Returns the <em>magnitude</em> of the angle between the first 
	 * and second blocks of {@code afpChain}, measured in degrees. 
	 * This is always a positive value (unsigned).
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	public static double getAngle(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {
		Matrix rotation = afpChain.getBlockRotationMatrix()[0];
		return Math.acos(rotation.trace() - 1) * 180/Math.PI;
	}

	/**
	 * Calculates a graph in the format of adjacency list from the 
	 * set of alignments, where each vertex is a residue and each 
	 * edge means the connection between the two residues in one 
	 * of the alignments.
	 * List dimensions: AdjList[vertices][edges]
	 * 
	 * @param allAFPs List of AFPChain
	 * @param ca1 Atom array of the symmetric structure
	 * @param undirected boolean make the graph undirected if true, 
	 * 			make it directed otherwise
	 * @return List double adjacency list defining a graph
	 */
	public static List<List<Integer>> buildAFPgraph(
			List<AFPChain> allAFPs, Atom[] ca1, 
			boolean undirected) {

		//Initialize the adjacency list that stores the graph
		List<List<Integer>> adjList = new ArrayList<List<Integer>>();
		for (int n=0; n<ca1.length; n++){
			List<Integer> edges = new ArrayList<Integer>();
			adjList.add(edges);
		}

		for (int k=0; k < allAFPs.size(); k++){
			for (int i=0; i<allAFPs.get(k).getOptAln().length; i++){
				for (int j=0; j<allAFPs.get(k).getOptAln()[i][0].length; j++){

					//The vertex is the residue in the first chain
					int vertex = allAFPs.get(k).getOptAln()[i][0][j];
					//The edge the one in the second chain
					int edge = allAFPs.get(k).getOptAln()[i][1][j];
					if (!adjList.get(vertex).contains(edge)){
						adjList.get(vertex).add(edge);
					}
					//Make the graph undirected (optional feature)
					if (undirected) {
						if (!adjList.get(edge).contains(vertex)) 
							adjList.get(edge).add(vertex);
					}
				}
			}
		}
		//Sort the edges in the adjacency list
		for (List<Integer> v:adjList) Collections.sort(v);
		return adjList;
	}

	/**
	 * Method that converts the symmetric units of a structure into different
	 * chains, so that internal symmetry can be translated into quaternary.<p>
	 * Application: obtain the internal symmetry axis with the quaternary 
	 * symmetry code in biojava.
	 * 
	 * @param symmetry MultipleAlignment of the subunits only
	 * @return Structure with different chains for every symmetric unit
	 */
	public static Structure toQuaternary(MultipleAlignment symmetry) {

		if (!symmetry.getEnsemble().getAlgorithmName().contains("symm")){
			throw new IllegalArgumentException(
					"The input alignment is not a symmetry alignment.");
		}

		Atom[] atoms = symmetry.getEnsemble().getAtomArrays().get(0);
		Structure cloned = atoms[0].getGroup().getChain().getParent().clone();
		atoms = StructureTools.getRepresentativeAtomArray(cloned);

		Structure symm = new StructureImpl();
		symm.setChains(new ArrayList<Chain>());
		char chainID = 'A';

		//Create new structure containing the subunit atoms
		for (int i=0; i<symmetry.size(); i++){
			Chain newCh = new ChainImpl();
			newCh.setChainID(chainID + "");
			chainID++;
			
			symm.addChain(newCh);
			Block align = symmetry.getBlocks().get(0);

			//Determine start and end of the subunit
			int count = 0;
			Integer start = null;
			while (start == null && count<align.length()){
				start = align.getAlignRes().get(i).get(0+count);
				count++;
			}
			count = 1;
			Integer end = null;
			while (end == null && count<symmetry.length()){
				end = align.getAlignRes().get(i).get(align.length()-count);
				count++;
			}
			end++;

			Atom[] subunit = Arrays.copyOfRange(atoms, start,end);

			for (int k=0; k<subunit.length; k++)
				newCh.addGroup((Group) subunit[k].getGroup().clone());
		}
		return symm;
	}

	/**
	 * Method that converts a subunit symmetric alignment into an alignment
	 * of whole structures.<p>
	 * Example: if the structure has subunits A,B and C, the original alignment
	 * is A-B-C, and the returned alignment is ABC-BCA-CAB.
	 * 
	 * @param symmetry MultipleAlignment of the subunits only
	 * @return Structure with different chains for every symmetric unit
	 */
	public static MultipleAlignment toFullAlignment(MultipleAlignment symm) {

		if (!symm.getEnsemble().getAlgorithmName().contains("symm")){
			throw new IllegalArgumentException(
					"The input alignment is not a symmetry alignment.");
		}

		MultipleAlignment full = symm.clone();

		for (int str=1; str<full.size(); str++){
			//Create a new Block with swapped AlignRes (move first to last)
			Block b = full.getBlocks().get(full.getBlocks().size()-1).clone();
			b.getAlignRes().add(b.getAlignRes().get(0));
			b.getAlignRes().remove(0);
			full.getBlockSets().get(0).getBlocks().add(b);
		}
		return full;
	}
}
