package org.biojava.nbio.structure.align.symm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Recursive version of CeSymm that aims at identifying all symmetry axis 
 * (internal or quaternary) of a particular structure.
 * <p>
 * Works in the following way:
 * <ul><li>Run CeSymm on the original structure.
 * <li>Calculate the symmetric unit boundaries.
 * <li>Run CeSymm on each of the symmetric units to find further symmetries.
 * <li>Repeat the last two steps until no more significant results are found.
 * <li>Run a final optimization of all symmetric units correctly superimposed.
 * </ul></li>
 * 
 * @author Aleix Lafita
 *
 */
public class CESymmRecursive {
	
	private CESymmParameters params;
	private MultipleAlignment msa;
	private int subunitLen;
	private Atom[] allAtoms;
	
	/**
	 * For the recursive algorithm to work properly the refinement and 
	 * optimization options should be turned on, because the alignment
	 * has to be consistent.
	 * 
	 * @param params CeSymm parameters
	 */
	public CESymmRecursive(CESymmParameters params){
		this.params = params;
		msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(new ArrayList<Atom[]>());
		msa.setTransformations(new ArrayList<Matrix4d>());
		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(new ArrayList<List<Integer>>());
		subunitLen = 0;
	}
	
	/**
	 * Uses a depth first search mode to calculate all symmetries in the input 
	 * array of atoms.
	 * After the calculation use {@link #getMultipleAlignment()} on the 
	 * recurser instance to get the result.
	 * 
	 * @param atoms
	 * @return used for the recursion only, ignore the return Integer
	 * @throws StructureException
	 */
	public int recurse(Atom[] atoms) throws StructureException {
		
		if (allAtoms == null) allAtoms = atoms;
		CeSymm aligner = new CeSymm();
		List<Atom[]> array = new ArrayList<Atom[]>();
		array.add(atoms);
		MultipleAlignment align = aligner.align(array, params);
		
		//Base case of the recursion
		if (align == null) return 0;
		else if (align.getScore(MultipleAlignmentScorer.AVGTM_SCORE) < 
				CeSymm.symmetryThreshold || align.length() < 20) {
			return 0;
		}
		
		Matrix4d t = align.getTransformations().get(align.size()-1);
		int children = 0;
		
		//If the alignment is meaningful, create one more level of recursion
		for (int su=0; su<align.size();su++) {
			
			Atom[] atomsR = Arrays.copyOfRange(atoms, 
					align.getBlocks().get(0).getAlignRes().get(su).get(0),
					align.getBlocks().get(0).getAlignRes().get(su).
					get(align.getBlocks().get(0).length()-1)+1);
						
			children = recurse(atomsR);
			
			if (children == 0){
				//This is the lowest level of recursion
				msa.getEnsemble().getAtomArrays().add(allAtoms);
				Matrix4d transform = new Matrix4d();
				transform.setIdentity();
				for (int b=0; b<su; b++) transform.mul(t);
				msa.getTransformations().add(transform);
				List<Integer> residues = 
						align.getBlocks().get(0).getAlignRes().get(su);
				msa.getBlocks().get(0).getAlignRes().add(residues);
				
				//Ensure that the size of the alignment is equal for all
				if (residues.size() > subunitLen) 
					subunitLen = residues.size();
				
				for (int i=0; i<msa.getBlocks().get(0).size(); i++){
					while (msa.getBlocks().get(0).getAlignRes().get(i).size() <
							subunitLen){
						msa.getBlocks().get(0).getAlignRes().get(i).add(null);
					}
				}
			} else {
				//If the lower levels succeeded just update the info
				Matrix4d transform = new Matrix4d();
				transform.setIdentity();
				for (int b=0; b<su; b++) {
					transform.mul(t);
				}
				for (int c = 0; c<children; c++){
					int index = msa.getTransformations().size() - children + c;
					msa.getTransformations().get(index).mul(transform);
					Block b = msa.getBlocks().get(0);
					int size = b.getAlignRes().get(index).size();
					for (int r=0; r<size; r++){
						Integer residue = b.getAlignRes().get(index).get(r);
						if (residue != null) {
							int sum = align.getBlocks().get(0).getAlignRes().
									get(su).get(0);
							b.getAlignRes().get(index).set(r, residue+sum);
						}
					}
				} //End of children update
			}
		}
		return align.size()+children;
	}
	
	/**
	 * Return the resulting symmetry alignment.
	 * @return MultipleAlignment
	 */
	public MultipleAlignment getMultipleAlignment(){
		return msa;
	}
	
	public static void main(String[] args) throws Exception {
		
		//Structures with more than one symmetry axis: 4gcr, 1vym.A, 1yox.A
		//Domain swapping: 1g6s
		//Internal+quaternary: 1VYM, 1f9z, 1YOX_A:,B:,C:
		//Structures that have different symmetry thresholds: 1vzw
		//Dihedral structures: 4hhb, 1iy9, 2ehz
		String name = "4gcr";
		
		AtomCache cache = new AtomCache();
		Atom[] atoms = ChainSorter.cyclicSorter(cache.getStructure(name));
		
		CESymmParameters params = new CESymmParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setOptimization(true);
		
		CESymmRecursive recurser = new CESymmRecursive(params);
		recurser.recurse(atoms);
		
		MultipleAlignment msa = recurser.getMultipleAlignment();
		
		msa.getEnsemble().setStructureNames(new ArrayList<String>());
		for (int s=0; s<msa.getBlocks().get(0).size(); s++) 
			msa.getEnsemble().getStructureNames().add(name+"_"+(s+1));
		
		new SymmetryJmol(msa);
	}
}
