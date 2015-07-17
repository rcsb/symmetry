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
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.axis.SymmetryAxes;
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
 * TODO waiting for jgrapht dependency to implement this
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmRecursive {

	private CESymmParameters params;
	private MultipleAlignment msa;
	private int subunitLen;
	private Atom[] allAtoms;
	private SymmetryAxes axes;
	//private Graph<MultipleAlignment> graph; TODO

	/**
	 * For the recursive algorithm to work properly the refinement and 
	 * optimization options should be turned on, because the alignment
	 * has to be consistent at every recursive step.
	 * 
	 * @param params CeSymm parameters
	 */
	public CeSymmRecursive(CESymmParameters params){
		this.params = params;

		msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(new ArrayList<Atom[]>());
		msa.getEnsemble().setAlgorithmName(CeSymm.algorithmName);
		msa.getEnsemble().setVersion(CeSymm.version);
		msa.getEnsemble().setStructureNames(new ArrayList<String>());

		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(new ArrayList<List<Integer>>());

		/*graph = new DirectedGraph<MultipleAlignment>();
		graph.addVertex(null); //root vertex*/

		axes = new SymmetryAxes();
		subunitLen = 0;
	}

	/**
	 * This method uses a DFS recursion mode to calculate all 
	 * symmetries in the input array of atoms and organize them in a
	 * connected graph of alignments.<p>
	 * Afterwards, it explores all vertices of the graph to determine
	 * the final MultipleAlignment of the subunits, taking care that 
	 * the recursion depth is the same for all subunits.
	 * 
	 * @param atoms atoms 
	 * @return MultipleAlignment of the subunits
	 * @throws Exception 
	 */
	public MultipleAlignment execute(Atom[] atoms) throws Exception {

		allAtoms = atoms;
		recurse(atoms, null);

		buildAlignment();
		return msa;
	}

	/**
	 * This method recursively creates the full Graph of symmetry
	 * alignments of the structure. The graph has to be analyzed 
	 * afterwards to determine the symmetric subunits.
	 * 
	 * @param atoms Coordinates of the symmetric structure
	 * @param parent parent vertex of the Graph
	 * @return int number of subunits found
	 * @throws StructureException
	 */
	private int recurse(Atom[] atoms, MultipleAlignment parent) 
			throws StructureException {

		//Perform the CeSymm alignment
		CeSymm aligner = new CeSymm();
		List<Atom[]> array = new ArrayList<Atom[]>();
		array.add(atoms);
		MultipleAlignment align = aligner.align(array, params);

		//Base case of the recursion, when the alignment is asymmetric
		/*if (align == null) return 1;
		else if (align.getScore(MultipleAlignmentScorer.AVGTM_SCORE) < 
				CeSymm.symmetryThreshold || align.length() < 20) {
			return 1;
		}

		//Case of symmetry found
		graph.addVertex(align);
		graph.addEdge(parent, align);

		int children = 1;

		//Create one more level of recursion
		for (int su=0; su<align.size(); su++) {

			//Generate the Atoms of the symmetric subunit
			Integer start = null;
			int it = 0;
			while (start == null){
				start = align.getBlocks().get(0).getAlignRes().get(su).get(it);
				it++;
			}
			Integer end = null;
			it = align.getBlocks().get(0).getAlignRes().get(su).size()-1;
			while (end == null){
				end = align.getBlocks().get(0).getAlignRes().get(su).get(it);
				it--;
			}

			Atom[] atomsR = Arrays.copyOfRange(atoms, start, end+1);
			children = recurse(atomsR, align);
		}

		return children*align.size();*/
		return 1;
	}

	private void buildAlignment() throws StructureException {

		//Find the depth in the graph that has equal levels of recursion
		/*int depth = 0;
		boolean stop = false;
		List<Integer> childNr = new ArrayList<Integer>();
		List<Integer> children = new ArrayList<Integer>();
		List<Integer> parents = new ArrayList<Integer>();
		children.add(0);
		while (!stop){
			childNr = new ArrayList<Integer>();
			parents = children;
			children = new ArrayList<Integer>();
			for (int p=0; p<parents.size(); p++){
				List<Integer> ch = graph.getChildren(parents.get(p));
				children.addAll(ch);
				childNr.add(ch.size());
			}
			//Check that no children is empty
			if (childNr.get(0) == 0) stop = true;
			else {
				//Check that all have the same number of children
				for (int i=1; i<childNr.size(); i++){
					if (childNr.get(0) != childNr.get(i)){
						stop = true;
						break;
					}
				}
			}
			if (childNr.size()==0) stop = true;
			depth++;
		}

		for (int index : parents){

			MultipleAlignment align = graph.getVertex(index);

			for (int su=0; su<align.size(); su++) {

				msa.getEnsemble().getAtomArrays().add(allAtoms);
				String name = align.getEnsemble().getStructureNames().get(0);
				msa.getEnsemble().getStructureNames().add(name);

				List<Integer> residues = 
						align.getBlocks().get(0).getAlignRes().get(su);

				//Sum the parent start residue recursively until the top
				int p = index;
				while (true){

					Integer c = p;
					p = graph.getParents(c).get(0);
					if (p == 0) break; //root reached

					MultipleAlignment parent = graph.getVertex(p);
					List<Integer> childrenOfParent = graph.getChildren(p);
					int subunit = childrenOfParent.indexOf(c);

					Integer start = null;
					int it = 0;
					while (start == null){
						start = parent.getBlocks().get(0).getAlignRes().
								get(subunit).get(it);
						it++;
					}
					for (int i=0; i<residues.size(); i++){
						Integer residue = residues.get(i);
						if (residue != null) residue += start;
						residues.set(i, residue);
					}
				}
				msa.getBlocks().get(0).getAlignRes().add(residues);

				//Ensure that the size of the alignment is equal for all
				if (residues.size() > subunitLen) {
					subunitLen = residues.size();
				}
				for (int i=0; i<msa.getBlocks().get(0).size(); i++){
					while (msa.getBlocks().get(0).getAlignRes().get(i).size() 
							< subunitLen){
						msa.getBlocks().get(0).getAlignRes().get(i).add(null);
					}
				}
			}
		}*/
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

		CeSymmRecursive recurser = new CeSymmRecursive(params);
		MultipleAlignment msa = recurser.execute(atoms);

		new SymmetryJmol(msa, recurser.getSymmetryAxes());
	}
}
