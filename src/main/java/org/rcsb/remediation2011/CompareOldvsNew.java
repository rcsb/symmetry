package org.rcsb.remediation2011;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.PDBStatus;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.io.LocalCacheStructureProvider;
import org.biojava.bio.structure.io.SandboxStyleStructureProvider;

public class CompareOldvsNew {

	List<String> nrAtomsChanged;
	List<String> caMissing;
	List<String> notTheSame;
	
	public static void main(String[] args){
		SandboxStyleStructureProvider sandbox2011 = new SandboxStyleStructureProvider();

		sandbox2011.setPath("/Users/ap3/WORK/PDB2011/");

		LocalCacheStructureProvider pdb2010 = new LocalCacheStructureProvider();
		pdb2010.setPath("/Users/ap3/WORK/PDB/");
		
		
		CompareOldvsNew compare = new CompareOldvsNew();

		try {
			List<String> currentPdbIds = sandbox2011.getAllPDBIDs();

			System.out.println("found " + currentPdbIds.size() + " PDB IDs");

			for ( String pdbId : currentPdbIds){

				compare.compare(pdbId,pdb2010,sandbox2011);


			}
		} catch (Exception e){
			e.printStackTrace();
		}

		
		System.out.println("The problem cases:");
		compare.printProblemCases();
	}
	
	private void printProblemCases() {
		System.out.println("Nr CA Atoms changed: " + nrAtomsChanged.size() + " " + nrAtomsChanged);
		System.out.println("CA atoms missing: " + caMissing.size() + " " + caMissing);
		System.out.println("Not the same CA, but same PDBresnum: " + notTheSame.size() + " " + notTheSame);
	}

	public CompareOldvsNew(){
		nrAtomsChanged = new ArrayList<String>();
		caMissing = new ArrayList<String>();
		notTheSame = new ArrayList<String>();
	}

	private void compare(String pdbId, LocalCacheStructureProvider pdb2010,
			SandboxStyleStructureProvider sandbox2011) throws IOException, StructureException {

		System.out.println("comparing: " + pdbId);

		PDBStatus.Status status = PDBStatus.getStatus(pdbId);
		if ( status == null)
			status = PDBStatus.Status.UNKNOWN;
		if ( ! status.equals(PDBStatus.Status.CURRENT)) {
			System.out.println("   not current...");
			return;
		}
		
		
		Structure oldS = pdb2010.getStructureById(pdbId);
		Structure newS = sandbox2011.getStructureById(pdbId);

		Atom[] newA = StructureTools.getAtomCAArray(newS);
		Atom[] oldA = StructureTools.getAtomCAArray(oldS);

		if ( newA.length != oldA.length){
			System.err.println("Length of CA atoms has changed: " + pdbId + " new: " + newA.length + " " + oldA.length);
			nrAtomsChanged.add(pdbId);
		}

		
		for (Chain newC : newS.getChains()){

			String chainId = newC.getChainID();

			if ( ! oldS.hasChain(chainId)) {
				System.err.println("Missing chain " + chainId + " in PDB " + pdbId);

				return;
			}

			Chain oldC = oldS.getChainByPDB(chainId);
			newA = StructureTools.getAtomCAArray(newC);
			oldA = StructureTools.getAtomCAArray(oldC);

			if ( newA.length != oldA.length){
				System.err.println("Length of CA atoms has changed: " + pdbId + " chain: " + chainId + "  new: " + newA.length + " " + oldA.length);
			}


			for ( int i = 0 ; i < oldA.length; i++){
				Atom a = oldA[i];
				Group oldG = a.getGroup();

				try {

					Group newG = newC.getGroupByPDB(oldG.getResidueNumber());

					Atom b = newG.getAtom(" CA ");
					double dist = Calc.getDistance(a, b);
					if ( dist > 0.001){
						System.err.println("Same PDBresnum, but not the same CA! " + pdbId + " " + oldG.getResidueNumber() + " " + dist);
						if ( ! notTheSame.contains(pdbId))
							notTheSame.add(pdbId);
					}
					
				} catch (Exception e){
					System.err.println(e.getMessage());
					if ( ! caMissing.contains(pdbId))
						caMissing.add(pdbId);
				}
			}

		}
	}

}
