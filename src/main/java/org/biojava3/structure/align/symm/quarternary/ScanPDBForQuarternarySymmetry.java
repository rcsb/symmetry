package org.biojava3.structure.align.symm.quarternary;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SimpleSubstitutionMatrix;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.util.ConcurrencyTools;
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.rcsb.fatcat.server.PdbChainKey;

public class ScanPDBForQuarternarySymmetry {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		AtomCache cache = new AtomCache();
		cache.setAutoFetch(true);

		int total = 0;
		Set<String> reps = getRepresentativePDBIds();
		System.out.println("Total representative PDBs: " + reps.size());
		
		long t1 = System.nanoTime();
		for (String pdbId : reps){
			System.out.println(pdbId);

			Structure structure = null;

			try {
				structure = cache.getBiologicalUnit(pdbId);
				structure.setPDBCode(pdbId); // pdb ids not set by default
			} catch (Exception e) {
				try {
					structure = cache.getStructure(pdbId);
				} catch (Exception e1) {
					System.err.println(e1.getMessage());
					continue;
				} 
			} 

			FindQuarternarySymmetry finder = new FindQuarternarySymmetry();
			finder.run(structure);
			
			total++;
			if (total == 100) break;
		}
	    long t2 = System.nanoTime();
	    System.out.println("Time: " + (t2-t1)/1000000 + " ms.");
	}

	private static Set<String> getRepresentativePDBIds() {
		System.out.println("GetRepresentatives");
		Set<String> set = new TreeSet<String>();
		SortedSet<PdbChainKey> reps = GetRepresentatives.getRepresentatives();
		for ( PdbChainKey r : reps){
			set.add(r.getPdbId());
		}
		return set;
	}
}
