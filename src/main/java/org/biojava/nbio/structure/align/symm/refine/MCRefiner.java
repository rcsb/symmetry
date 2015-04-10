package org.biojava.nbio.structure.align.symm.refine;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.subunit.MultipleAFP;
import org.biojava.nbio.structure.align.util.AlignmentTools;

/**
 * Creates a refined alignment by a MC optimization of the subunit multiple alignment.
 * @author lafita
 */

public class MCRefiner implements Refiner {

	public MCRefiner() {
		super();
	}
	
	@Override
	public AFPChain refine(AFPChain[] afpAlignments, Atom[] ca1, Atom[] ca2, int order)
			throws RefinerFailedException,StructureException {
		
		AFPChain originalAFP = afpAlignments[0];
		AFPChain refinedAFP = SingleRefiner.refineSymmetry(originalAFP, ca1, ca2, order);
		
		AFPChain finalAFP = partitionAFPchain(refinedAFP,ca1,ca2,order);
			
		MultipleAFP mulAln = new MultipleAFP(finalAFP,ca1);
		
		//return finalAFP;
		return mulAln.getAfpChain();
	}
	
	/**
	 *  Partitions an afpChain alignment into order blocks of aligned residues of the same size.
	 */
	private AFPChain partitionAFPchain(AFPChain afpChain, Atom[] ca1, Atom[] ca2, int order) throws StructureException{
		
		int[][][] newAlgn = new int[order][][];
		int subunitLen = afpChain.getOptLength()/order;
		
		//Extract all the residues considered in the first chain of the alignment
		List<Integer> alignedRes = new ArrayList<Integer>();
		for (int su=0; su<afpChain.getBlockNum(); su++){
			for (int i=0; i<afpChain.getOptLen()[su]; i++){
				alignedRes.add(afpChain.getOptAln()[su][0][i]);
			}
		}
		
		//Build the new alignment, leaving two residue spacing between subunits to avoid hard boundaries 
		//and two residue spacing in the middle of the subunits to allow shifts during optimization
		for (int su=0; su<order; su++){
			newAlgn[su] = new int[2][];
			newAlgn[su][0] = new  int[subunitLen];
			newAlgn[su][1] = new  int[subunitLen];
			for (int i=0; i<subunitLen; i++){
				newAlgn[su][0][i] = alignedRes.get(subunitLen*su+i);
				newAlgn[su][1][i] = alignedRes.get((subunitLen*(su+1)+i)%alignedRes.size());
			}
		}
		
		return AlignmentTools.replaceOptAln(newAlgn, afpChain, ca1, ca2);
	}
}