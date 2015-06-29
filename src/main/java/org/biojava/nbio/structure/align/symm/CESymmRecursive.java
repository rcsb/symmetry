package org.biojava.nbio.structure.align.symm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.utils.SymmetryTools;

/**
 * Recursive version of CeSymm that aims at identifying all symmetry axis (internal or quaternary) 
 * of a particular structure.
 * <p>
 * Works in the following way:
 * <ul><li>Run CeSymm on the original structure.
 * <li>Calculate the symmetric unit boundaries.
 * <li>Run CeSymm on each of the symmetric units to find further symmetries.
 * <li>Repeat the last two steps until no more significant results are found.
 * <li>Combine all the symmetry axis to find the smallest subset of axis that describe the symmetry of the structure.
 * <li>Run a final optimization of all symmetric units correctly superimposed.
 * </ul></li>
 * The recursion levels are calculated in different threads in order to improve the running time.
 * 
 * 
 * @author Aleix Lafita
 *
 */
public class CESymmRecursive {
	
	private CESymmParameters params;
	private MultipleAlignment msa;
	private List<Matrix4d> axis;
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
		axis = new ArrayList<Matrix4d>();
		subunitLen = 0;
	}
	
	/**
	 * Uses a depth first search mode to calculate all symmetries in the input array of atoms.
	 * After the calculation use {@link #getMultipleAlignment()} on the recurser instance
	 * to get the result.
	 * 
	 * @param atoms
	 * @return used for the recursion only, ignore the return Integer
	 * @throws StructureException
	 */
	public int recurse(Atom[] atoms) throws StructureException {
		
		if (allAtoms == null) allAtoms = atoms;
		CeSymm aligner = new CeSymm();
		AFPChain afp = aligner.align(atoms, atoms, params);
		
		//Base case of the recursion
		if (afp.getTMScore() < CeSymm.symmetryThreshold || afp.getOptLength() < 8) return 0;
		
		//Calculate the symmetry axis and add it to the list only if it is unique
		Matrix4d matrixT = Calc.getTransformation(afp.getBlockRotationMatrix()[0], afp.getBlockShiftVector()[0]);
		boolean equivalent = false;
		for (int i=0; i<axis.size(); i++){
			if (!SymmetryTools.areEquivalentAxis(matrixT, axis.get(i), 0.2)) continue;
			else {
				equivalent = true;
				break;
			}
		}
		if (!equivalent) axis.add(matrixT);
		
		//If the alignment is meaningful, create one more level of recursion
		for (int bk=0; bk<afp.getBlockNum(); bk++){
			Atom[] atomsR = Arrays.copyOfRange(atoms, afp.getOptAln()[bk][0][0], afp.getOptAln()[bk][0][afp.getOptLen()[bk]-1]+1);
						
			int children = recurse(atomsR);
			if (children == 0){
				//This is the lowest level of recursion, so add the residues and matrix information
				msa.getEnsemble().getAtomArrays().add(allAtoms);
				Matrix4d transform = new Matrix4d();
				transform.setIdentity();
				for (int b=0; b<bk; b++) transform.mul(matrixT);
				msa.getTransformations().add(transform);
				List<Integer> residues = new ArrayList<Integer>();
				for (Integer res:afp.getOptAln()[bk][0]) residues.add(res);
				msa.getBlocks().get(0).getAlignRes().add(residues);
				
				//Ensure that the size of the alignment is equal for all subunits. Insert gaps otherwise
				if (residues.size() > subunitLen) subunitLen = residues.size();
				for (int i=0; i<msa.getBlocks().get(0).size(); i++){
					while (msa.getBlocks().get(0).getAlignRes().get(i).size() < subunitLen){
						msa.getBlocks().get(0).getAlignRes().get(i).add(null);
					}
				}
			}
			else {
				//If the lower levels succeeded just update the matrices of the childs
				Matrix4d transform = new Matrix4d();
				transform.setIdentity();
				for (int b=0; b<bk; b++) transform.mul(matrixT);
				for (int child = 0; child<children; child++){
					int index = msa.getTransformations().size() - children + child;
					msa.getTransformations().get(index).mul(transform);
					for (int res=0; res<msa.getBlocks().get(0).getAlignRes().get(index).size(); res++){
						Integer residue = msa.getBlocks().get(0).getAlignRes().get(index).get(res);
						if (residue != null) {
							msa.getBlocks().get(0).getAlignRes().get(index).set(res, residue+afp.getOptAln()[bk][0][0]);
						}
					}
				}
			}
		}
		return afp.getBlockNum();
	}
	
	/**
	 * Return the resulting symmetry alignment.
	 * @return MultipleAlignment
	 */
	public MultipleAlignment getMultipleAlignment(){
		return msa;
	}
	
	/**
	 * Return all the symmetry axis found in the structure.
	 * @return List of RotationAxis
	 */
	public List<Matrix4d> getAxis(){
		return axis;
	}
	
	public static void main(String[] args) throws Exception {
		
		//Structures with more than one symmetry axis: 4gcr, 1vym.A
		//Same symmetry axis but combined internal/quaternary: 1VYM, 1F9Z
		//Structures that have different symmetry thresholds: 1vzw
		//Dihedral structures: 4hhb, 1iy9
		String name = "1VYM.A";
		
		AtomCache cache = new AtomCache();
		Atom[] atoms = ChainSorter.cyclicSorter(cache.getStructure(name));
		
		CESymmParameters params = new CESymmParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setOptimization(true);
		
		CESymmRecursive recurser = new CESymmRecursive(params);
		recurser.recurse(atoms);
		
		MultipleAlignment msa = recurser.getMultipleAlignment();
		msa.getEnsemble().setStructureNames(new ArrayList<String>());
		for (int s=0; s<msa.getBlocks().get(0).size(); s++) msa.getEnsemble().getStructureNames().add(name+"_"+(s+1));
		
		//MultipleAlignmentOptimizerMC optimizer = new MultipleAlignmentOptimizerMC(msa, new MultipleMcParameters(), 0);
		//msa = optimizer.call();
		
		List<RotationAxis> symmAxis = new ArrayList<RotationAxis>();
		for (int i=0; i<recurser.getAxis().size(); i++){
			RotationAxis axis = new RotationAxis(recurser.getAxis().get(i));
			symmAxis.add(axis);
		}
		new SymmetryJmol(msa,symmAxis);
	}
}
