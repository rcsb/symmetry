package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
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
import org.biojava3.structure.align.symm.protodomain.Protodomain;
import org.biojava3.structure.align.symm.protodomain.ProtodomainCreationException;

/**
 * Find sequence motifs that correspond to symmetry subunits or protodomains.
 * @author dmyersturnbull
 */
public class ConservedSequenceFinder {

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2) {
			System.err.println("Usage: " + ConservedSequenceFinder.class.getSimpleName() + " input-census-file.xml [min-score]");
			return;
		}
		float minScore = 0.8f;
		if (args.length > 1) {
			minScore = Float.parseFloat(args[1]);
		}
		ConservedSequenceFinder finder = new ConservedSequenceFinder(new AtomCache());
		finder.find(Results.fromXML(new File(args[0])), minScore);
	}
	
	private AtomCache cache;
	
	public ConservedSequenceFinder(AtomCache cache) {
		super();
		this.cache = cache;
	}

	public void find(Results results, float maxDistance) {
		for (Result result : results.getData()) {
			if (result.getOrder() != 2) continue;
			Protodomain wholeAligned, protodomain1, protodomain2;
			try {
				wholeAligned = Protodomain.fromString(result.getProtodomain(), result.getScopId(), cache);
				protodomain1 = wholeAligned.createSubstruct(result.getOrder(), 0).spliceApproxConsecutive();
				protodomain2 = wholeAligned.createSubstruct(result.getOrder(), 1).spliceApproxConsecutive();
				System.out.println(wholeAligned);
				System.out.println(protodomain1);
				System.out.println(protodomain2);
				System.out.println();
			} catch (ProtodomainCreationException e) {
				e.printStackTrace();
				continue;
			}
			try {
				protodomain1.buildStructure();
				protodomain2.buildStructure();
			} catch (IOException e) {
				e.printStackTrace();
				continue;
			} catch (StructureException e) {
				e.printStackTrace();
				continue;
			}
			Atom[] ca1 = StructureTools.getAtomCAArray(protodomain1.getStructure());
			Atom[] ca2 = StructureTools.getAtomCAArray(protodomain2.getStructure());
			ProteinSequence seq1 = new ProteinSequence(StructureTools.convertAtomsToSeq(ca1));
			ProteinSequence seq2 = new ProteinSequence(StructureTools.convertAtomsToSeq(ca2));
			PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments.getPairwiseAligner(seq1, seq2, PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(), SubstitutionMatrixHelper.getBlosum62());
			SequencePair<ProteinSequence, AminoAcidCompound> pair = aligner.getPair();
			if (aligner.getDistance() < maxDistance) continue;
			// TODO we are really only interested in motifs in the equivalent residues in the protodomains
			System.out.println(pair.toString(80));
		}
	}
	
}
