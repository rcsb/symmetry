package org.biojava3.structure.align.symm.census2.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

/**
 * Find sequence motifs that correspond to symmetry subunits or protodomains.
 * @author dmyersturnbull
 */
public class ConservedSequenceFinder {

	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 3) {
			System.err.println("Usage: " + ConservedSequenceFinder.class.getSimpleName() + " input-census-file.xml [min-score]");
			return;
		}
		float minScore = 0.7f;
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(args[1]))));
		if (args.length > 2) {
			minScore = Float.parseFloat(args[1]);
		}
		ConservedSequenceFinder finder = new ConservedSequenceFinder(new AtomCache());
		finder.find(Results.fromXML(new File(args[0])), minScore, pw);
		pw.close();
	}

	private AtomCache cache;

	public ConservedSequenceFinder(AtomCache cache) {
		super();
		this.cache = cache;
	}

	public void find(Results results, float maxDistance, PrintWriter pw) {
		
		Collections.shuffle(results.getData()); // shuffle!
		
		for (Result result : results.getData()) {
			
			if (!SignificanceFactory.rotationallySymmetricSmart().isSignificant(result)) continue;
			
			if (result.getOrder() == null || result.getOrder() != 2) continue;
			
			Protodomain wholeAligned, protodomain1, protodomain2;
			
			try {
				
				wholeAligned = Protodomain.fromString(result.getProtodomain(), result.getScopId(), cache);
				protodomain1 = wholeAligned.createSubstruct(result.getOrder(), 0).spliceApproxConsecutive();
				protodomain2 = wholeAligned.createSubstruct(result.getOrder(), 1).spliceApproxConsecutive();
				protodomain1.buildStructure();
				protodomain2.buildStructure();
				
				Atom[] ca1 = StructureTools.getAtomCAArray(protodomain1.getStructure());
				Atom[] ca2 = StructureTools.getAtomCAArray(protodomain2.getStructure());
				
				ProteinSequence seq1 = new ProteinSequence(StructureTools.convertAtomsToSeq(ca1));
				ProteinSequence seq2 = new ProteinSequence(StructureTools.convertAtomsToSeq(ca2));
				
				PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments.getPairwiseAligner(seq1, seq2, PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(), SubstitutionMatrixHelper.getBlosum62());
				SequencePair<ProteinSequence, AminoAcidCompound> pair = aligner.getPair();
//				System.err.println(aligner.getDistance());
				
				// TODO we are really only interested in motifs in the equivalent residues in the protodomains
				if (aligner.getDistance() > maxDistance) continue;
				
				pw.println(result.getScopId());
				pw.println(protodomain1);
				pw.println(protodomain2);
				
				pw.println(pair.toString(80));
				pw.println("score: " + aligner.getScore());

				pw.println("-----------------------------------------------------");
				
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
			}
		}
	}

}
