package org.biojava.nbio.structure.align.symm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.axis.SymmetryAxes;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.symmetry.utils.DirectedGraph;
import org.biojava.nbio.structure.symmetry.utils.Graph;

/**
 * Iterative version of CeSymm that aims at identifying all symmetry axis 
 * (internal or quaternary) of a particular structure.
 * <p>
 * Works in the following way:
 * <ul><li>Run CeSymm on the original structure.
 * <li>Calculate the symmetric unit boundaries.
 * <li>Run CeSymm on one of the symmetric units to find further symmetries.
 * <li>Repeat the last two steps until no more significant results are found.
 * <li>Map back all residues in a multiple alignment of the subunits.
 * <li>Run a final optimization of all symmetric units correctly superimposed.
 * </ul></li>
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmIterative {

	private CESymmParameters params;
	private MultipleAlignment msa;
	private Atom[] allAtoms;
	private String name;
	private SymmetryAxes axes;
	private Graph<Integer> alignment;

	/**
	 * For the iterative algorithm to work properly the refinement and 
	 * optimization options should be turned on, because the alignment
	 * has to be consistent at every recursive step.
	 * 
	 * @param params CeSymm parameters
	 */
	public CeSymmIterative(CESymmParameters params) {
		this.params = params;

		msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(new ArrayList<Atom[]>());
		msa.getEnsemble().setAlgorithmName(CeSymm.algorithmName);
		msa.getEnsemble().setVersion(CeSymm.version);
		msa.getEnsemble().setStructureNames(new ArrayList<String>());

		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(new ArrayList<List<Integer>>());

		alignment = new DirectedGraph<Integer>();
		axes = new SymmetryAxes();
		name = null;
	}

	/**
	 * This method uses iteratively CeSymm to calculate all 
	 * symmetries in the input array of atoms and organize 
	 * them in a multiple alignment of the subunits.
	 * 
	 * @param atoms atoms 
	 * @return MultipleAlignment of the subunits
	 * @throws Exception 
	 */
	public MultipleAlignment execute(Atom[] atoms) throws Exception {

		allAtoms = atoms;
		for (Integer res=0; res<allAtoms.length; res++){
			alignment.addVertex(res);
		}

		iterate(atoms, 0);
		buildAlignment();

		//Run a final optimization of the alignment
		/*
		MultipleMcParameters params = new MultipleMcParameters();
		params.setMinAlignedStructures(2);
		params.setMinBlockLen(15);
		MultipleMcOptimizer optimizer = 
				new MultipleMcOptimizer(msa, params, 0);

		//Not yet possible TODO 
		msa = optimizer.call();

		//Reset the atom Arrays for consistency with the symmetry format
		msa.getEnsemble().setAtomArrays(new ArrayList<Atom[]>());
		for (int su=0; su<msa.size(); su++){
			msa.getEnsemble().getAtomArrays().add(allAtoms);
		}*/
		return msa;
	}

	/**
	 * This method runs iteratively CeSymm on the symmetric units
	 * until no more symmetries exist.
	 * 
	 * @param atoms Coordinates of the symmetric structure
	 * @param first starting position of the atom array in the original array
	 * @throws StructureException
	 */
	private void iterate(Atom[] atoms, int first) throws StructureException {

		//Perform the CeSymm alignment
		CeSymm aligner = new CeSymm();
		List<Atom[]> array = new ArrayList<Atom[]>();
		array.add(atoms);
		MultipleAlignment align = aligner.align(array, params);
		if (name == null) 
			name = align.getEnsemble().getStructureNames().get(0);

		//End iterations if non symmetric
		if (align == null) return;
		else if (align.getScore(MultipleAlignmentScorer.AVGTM_SCORE) < 
				CeSymm.symmetryThreshold || align.length() < 20) {
			return;
		}

		//If symmetric store the residue dependencies in graph
		Block b = align.getBlocks().get(0);
		for (int pos=0; pos<b.length(); pos++){
			for (int su=0; su<b.size()-1; su++){
				Integer pos1 = b.getAlignRes().get(su).get(pos);
				Integer pos2 = b.getAlignRes().get(su+1).get(pos);
				//Add edge from lower to higher positions
				if (pos1 != null && pos2 != null){
					alignment.addEdge(pos1, pos2);
				}
			}
		}

		//Generate the Atoms of one of the symmetric subunit
		Integer start = null;
		int it = 0;
		while (start == null){
			start = align.getBlocks().get(0).getAlignRes().get(0).get(it);
			it++;
		}
		Integer end = null;
		it = align.getBlocks().get(0).getAlignRes().get(0).size()-1;
		while (end == null){
			end = align.getBlocks().get(0).getAlignRes().get(0).get(it);
			it--;
		}
		//Iterate further
		Atom[] atomsR = Arrays.copyOfRange(atoms, start, end+1);
		iterate(atomsR, start+first);
	}

	private void buildAlignment() throws StructureException {

		List<List<Integer>> groups = new ArrayList<List<Integer>>();
		List<Integer> alreadySeen = new ArrayList<Integer>();
		int size = 0;

		//Calculate the connected groups of the alignment graph
		for (int i=0; i<alignment.size(); i++){
			if (!alreadySeen.contains(i)){
				List<Integer> group = new ArrayList<Integer>();
				List<Integer> residues = new ArrayList<Integer>();
				residues.add(i);

				while (residues.size() > 0){
					List<Integer> newResidues = new ArrayList<Integer>();
					for (Integer residue : residues){
						group.add(residue);
						alreadySeen.add(residue);
						List<Integer> children = 
								alignment.getChildren(residue);
						newResidues.addAll(children);
					}
					residues = newResidues;
				}
				Collections.sort(group);
				groups.add(group);
				if (group.size() > size) size = group.size();
			}
		}
		
		Block b = msa.getBlocks().get(0);
		//Construct the MultipleAlignment
		for (int su=0; su<size; su++) {
			msa.getEnsemble().getStructureNames().add(name);
			msa.getEnsemble().getAtomArrays().add(allAtoms);
			b.getAlignRes().add(new ArrayList<Integer>());
			
			for (List<Integer> group : groups){
				if (group.size() != size) continue;
				b.getAlignRes().get(su).add(group.get(su));
			}
		}
	}

	/**
	 * Return the symmetry axes.
	 * @return SymmetryAxes
	 */
	public SymmetryAxes getSymmetryAxes(){
		return axes;
	}

	public static void main(String[] args) throws Exception {

		//More than one symmetry axis: 4gcr, 1ppr.O, 1vym.A, 1yox.A
		//Domain swapping: 1g6s
		//Internal+quaternary: 1VYM, 1f9z, 1YOX_A:,B:,C:, 1mmi
		//Structures that have different symmetry thresholds: 1vzw
		//Dihedral structures: 4hhb, 1iy9, 2ehz,
		String name = "1g6s";

		AtomCache cache = new AtomCache();
		Atom[] atoms = ChainSorter.cyclicSorter(cache.getStructure(name));

		CESymmParameters params = new CESymmParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setOptimization(true);

		CeSymmIterative aligner = new CeSymmIterative(params);
		MultipleAlignment msa = aligner.execute(atoms);

		new SymmetryJmol(msa, aligner.getSymmetryAxes());
	}
}
