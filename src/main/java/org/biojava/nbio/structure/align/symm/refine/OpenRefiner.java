package org.biojava.nbio.structure.align.symm.refine;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.gui.SymmetryJmol;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.utils.SymmetryTools;

/**
 * Creates an open symmetry refined alignment with the self-alignments obtained from CeSymm.
 * The goal is to find the internal repeats without the closed symmetry constraints.
 * The symmetry axis indicates also the translation of the subunits in this case.
 * It does not work for close symmetry, because it assumes that no CP ocurred in the first alignment.
 * 
 * @author Aleix Lafita
 * 
 */
public class OpenRefiner implements Refiner {

	public OpenRefiner() {
		super();
	}
	
	@Override
	public AFPChain refine(List<AFPChain> afpAlignments, Atom[] ca1, Atom[] ca2, int order) throws RefinerFailedException,StructureException {
		
		//Create a directed graph from the alignments
		List<List<Integer>> graph = SymmetryTools.buildAFPgraph(afpAlignments, ca1, true);
		AFPChain afpChain = afpAlignments.get(afpAlignments.size()-1);
		List<Integer> alreadySeen = new ArrayList<Integer>();
		
		//Calculate the connected groups of the alignment graph
		List<List<Integer>> groups = new ArrayList<List<Integer>>();
		for (int i=0; i<graph.size(); i++){
			if (!alreadySeen.contains(i)){
				List<Integer> group = new ArrayList<Integer>();
				int residue = i;
				while (residue != -1 && !alreadySeen.contains(residue)){
					group.add(residue);
					//Go to the next residue in sequence: two vertices mean (previous, next) in the graph
					if (graph.get(residue).size() > 1){
						if (graph.get(residue).get(1) > residue) residue = graph.get(residue).get(1);
						else residue = -1;
					}
					else if (graph.get(residue).size() > 0) {
						if (graph.get(residue).get(0) > residue) residue = graph.get(residue).get(0);
						else residue = -1;
					}
					else residue = -1;  //This means that the residue does not have a next
				}
				Collections.sort(group);
				if (group.size() > 1){
					groups.add(group);
					alreadySeen.addAll(group);
				}
			}
		}
		
		//Calculate the most common group size
		List<Integer> sizes = new ArrayList<Integer>(ca1.length);
		for (int p=0; p<ca1.length; p++) sizes.add(0);
		
		for (int i=0; i<groups.size(); i++){
			int gorder = groups.get(i).size();
			sizes.set(gorder, sizes.get(gorder)+1);
		}
		int maxNr = 0; //the total number of residues aligned of the subunits - max determines the order
		for (int s=2; s<sizes.size(); s++){
			if (sizes.get(s)*s > maxNr) {
				order = s;
				maxNr = sizes.get(s)*s;
			}
		}
		
		//Now create the new AFP alignment from the selected groups
		List<List<Integer>> subunits = new ArrayList<List<Integer>>();
		for (List<Integer> g:groups) if (g.size() == order) subunits.add(g); //add the groups with the right order
		
		//Delete all inconsistent groups in subunits (if they define different subunits)
		List<Integer> deleteIndices = new ArrayList<Integer>();
		for (int i=1; i<subunits.size(); i++){
			for (int j=0; j<subunits.get(i).size()-1; j++){
				if (subunits.get(i).get(j) > subunits.get(0).get(j+1)){
					deleteIndices.add(i);
					break;
				}
			}
		}
		for (int i=deleteIndices.size()-1; i>=0; i--){
			int index = deleteIndices.get(i);
			subunits.remove(index);
		}
		
		//From the groups of higher order take the consistent residues only (the ones that fall inside the subunits) - needs review
		for (List<Integer> g:groups){
			if (g.size() > order){
				List<Integer> group = new ArrayList<Integer>();
				for (int pos=0; pos<g.size() && group.size() < order; pos++){
					boolean consistent = true;
					for (List<Integer> sub:subunits){
						if (sub.get(group.size()) > g.get(pos)) consistent = false;
						if (group.size() < order-1) {
							if (sub.get(group.size()+1) < g.get(pos)) consistent = false;
						}
					}
					if (consistent && group.size()<order) group.add(g.get(pos));
				}
				if (group.size()==order) subunits.add(group);
			}
		}
		
		int[][][] optAln = new int[order][2][subunits.size()];
		for (int bk=0; bk<order; bk++){
			optAln[bk] = new int[2][];
			optAln[bk][0] = new int[subunits.size()];
			optAln[bk][1] = new int[subunits.size()];
			for (int su=0; su<subunits.size(); su++){
				optAln[bk][0][su] = subunits.get(su).get(bk);
				optAln[bk][1][su] = subunits.get(su).get((bk+1)%order);
			}
		}
		
		//Replace the alignment information without changing the superimposition
		afpChain = AlignmentTools.replaceOptAln(optAln, afpChain, ca1, ca2);
		
		return afpChain;
	}
	
	public static void main(String[] args) throws StructureException, IOException{
		
		String name = "1ppr.O";  //Ankyrin: 1N0R.A, 3EU9.A, 1AWC.B, 3EHQ.A, 1NFI.E
								  //Helical: 1EZG.A, 1D0B.A
								  //LRR: 2bnh.A, 1dfj.I
								  //HEAT: 1B3U.A
								  //TPR: 2FO7
								  //Ig-repeats: 2rik.A, 1fnf.A
								  //hevein: 1k7u.A
								  //Rossman fold: d1heta2
								  //benchmark H: d1af0a1, d1dcec3, d1kx9a_, d1l0sa_
								  //benchmark NIH: d1v0fd2, d2b1ea_
								  //benchmark SH: d1qtea1, d2ajab1
								  //benchmark R: d1blxb_, d1rmga_, d3bsda_
								  //benchmark C6: d1wp5a_
								  //TIM barrel duplication and twist: 1pii
								  //beta-hairpin: d2biba1
								  //duplication: 1vym.A, 1ppr.O
								  //clear repeats-turn-repeats example: 1S70.B - CeSymm makes slip alignment near identity...

		AtomCache cache = new AtomCache();

		//Parse atoms of the protein into two DS
		Atom[] ca1 = cache.getAtoms(name);
		Atom[] ca2 = cache.getAtoms(name);
		
		//Initialize a new CeSymm class and its parameters and a new alignment class
		CeSymm ceSymm = new CeSymm();
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setSymmetryType(SymmetryType.OPEN);
		params.setOptimization(true);
		//params.setSeed(11);
		AFPChain afpChain = new AFPChain();
		
		//Perform the alignment and store
		afpChain = ceSymm.align(ca1, ca2);
		afpChain.setName1(name);
		afpChain.setName2(name);
		
		//Display the AFP alignment of the subunits
		SymmetryJmol jmol = new SymmetryJmol(afpChain, ca1);
		//StructureAlignmentJmol jmol2 = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
		jmol.setTitle(name);
	}
}
