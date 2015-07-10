package org.biojava.nbio.structure.align.symm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcOptimizer;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcParameters;
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
		msa.getEnsemble().setAlgorithmName(CeSymm.algorithmName);
		msa.getEnsemble().setVersion(CeSymm.version);
		msa.getEnsemble().setStructureNames(new ArrayList<String>());

		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(new ArrayList<List<Integer>>());
		subunitLen = 0;
	}
	
	/**
	 * Uses a depth first search mode to calculate all symmetries in the input 
	 * array of atoms.
	 * 
	 * @param atoms
	 * @return used for the recursion only, ignore the return Integer
	 * @throws Exception 
	 */
	public MultipleAlignment execute(Atom[] atoms) throws Exception {
		
		allAtoms = atoms;
		recurse(atoms);
		
		//Run a final optimization of the alignment
		MultipleMcParameters params = new MultipleMcParameters();
		params.setMinAlignedStructures(2);
		params.setMinBlockLen(15);
		MultipleMcOptimizer optimizer = 
				new MultipleMcOptimizer(msa, params, 0);
		
		//Not yet possible TODO msa = optimizer.call();
		
		//Reset the atom Arrays for consistency with the symmetry format
		msa.getEnsemble().setAtomArrays(new ArrayList<Atom[]>());
		for (int su=0; su<msa.size(); su++){
			msa.getEnsemble().getAtomArrays().add(allAtoms);
		}
		return msa;
	}

	public int recurse(Atom[] atoms) throws StructureException {

		CeSymm aligner = new CeSymm();
		List<Atom[]> array = new ArrayList<Atom[]>();
		array.add(atoms);
		MultipleAlignment align = aligner.align(array, params);

		//Base case of the recursion
		if (align == null) return 1;
		else if (align.getScore(MultipleAlignmentScorer.AVGTM_SCORE) < 
				CeSymm.symmetryThreshold || align.length() < 20) {
			return 1;
		}

		int children = 1;
		int total = 0;

		//If the alignment is meaningful, create one more level of recursion
		for (int su=0; su<align.size(); su++) {

			int start = align.getBlocks().get(0).getAlignRes().get(su).get(0);
			int end = align.getBlocks().get(0).getAlignRes().get(su).get(
					align.getBlocks().get(0).getAlignRes().get(su).size()-1)+1;
			Atom[] atomsR = Arrays.copyOfRange(atoms, start, end);
			
			children = recurse(atomsR);
			total += children;

			if (children == 1){
				//This is the lowest level of recursion
				msa.getEnsemble().getAtomArrays().add(allAtoms);
				String name = align.getEnsemble().getStructureNames().get(0);
				msa.getEnsemble().getStructureNames().add(name);
				
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
				for (int c = 0; c<children; c++){
					int index = msa.size() - children + c;
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
		//TODO check that total == align.size()*children
		//or alternatively that nodes have the same children/recursion level
		return total;
	}

	/**
	 * Return the resulting symmetry alignment.
	 * @return MultipleAlignment
	 */
	public MultipleAlignment getMultipleAlignment(){
		return msa;
	}

	public static void main(String[] args) throws Exception {

		//More than one symmetry axis: 4gcr, 1ppr.O, 1vym.A, 1yox.A
		//Domain swapping: 1g6s
		//Internal+quaternary: 1VYM, 1f9z, 1YOX_A:,B:,C:
		//Structures that have different symmetry thresholds: 1vzw
		//Dihedral structures: 4hhb, 1iy9, 2ehz,
		String name = "1YOX_A:,B:,C";

		AtomCache cache = new AtomCache();
		Atom[] atoms = ChainSorter.cyclicSorter(cache.getStructure(name));

		CESymmParameters params = new CESymmParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setOptimization(true);

		CESymmRecursive recurser = new CESymmRecursive(params);
		MultipleAlignment msa = recurser.execute(atoms);

		new SymmetryJmol(msa);
	}
}
