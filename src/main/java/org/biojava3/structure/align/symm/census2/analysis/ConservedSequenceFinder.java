package org.biojava3.structure.align.symm.census2.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.SortedSet;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.SubstitutionMatrixScorer;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.protodomain.Protodomain;
import org.biojava3.ws.hmmer.HmmerResult;
import org.biojava3.ws.hmmer.HmmerScan;
import org.biojava3.ws.hmmer.RemoteHmmerScan;

/**
 * Find sequence motifs that correspond to symmetry subunits or protodomains.
 * @author dmyersturnbull
 */
public class ConservedSequenceFinder {

	private static final Logger logger = LogManager.getLogger(ConservedSequenceFinder.class.getName());

	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 3) {
			System.err.println("Usage: " + ConservedSequenceFinder.class.getSimpleName() + " input-census-file.xml [min-score]");
			return;
		}
		float minScore = 0.8f;
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(args[1]))));
		if (args.length > 2) {
			minScore = Float.parseFloat(args[1]);
		}
		ConservedSequenceFinder finder = new ConservedSequenceFinder(new AtomCache());
		finder.find(Results.fromXML(new File(args[0])), minScore, pw);
		pw.close();
	}

	private AtomCache cache;

	private double pValue;

	public ConservedSequenceFinder(AtomCache cache) {
		super();
		this.cache = cache;
	}

	private HmmerResult getPfam(ProteinSequence sequence) throws IOException {
		HmmerScan hmmer = new RemoteHmmerScan();
		SortedSet<HmmerResult> results = hmmer.scan(sequence);
		logger.debug("Found " + results.size() + " HMMER results");
		for (HmmerResult result : results) {
			logger.debug("Found result " + result);
			if (result.getPvalue() <= this.pValue) {
				logger.debug("Returning result " + result.getName());
				return result;
			}
		}
		return null;
	}

	/**
	 * Scores with a gamma distribution as per Webber and Barton 2001, Bioinformatics. The authors used (among other scoring
	 * schemes) the BLOSUM62 matrix with a gap opening penalty of 12 and extension penalty of 1. The paper is <a
	 * href="http://bioinformatics.oxfordjournals.org/content/17/12/1158.full.pdf+html">available</a>.
	 * 
	 * @author dmyersturnbull
	 * 
	 */
	public static class GammaScorer {
		private double alpha; // shape (normally theta)
		private double beta; // scale
		private double lambda;

		public static GammaScorer forBlosum62() {
			return new GammaScorer(25.54, 4.96, 0.2);
		}

		public GammaScorer(double alpha, double beta, double lambda) {
			super();
			this.alpha = alpha;
			this.beta = beta;
			this.lambda = lambda;
		}

		public double score(SequencePair<ProteinSequence, AminoAcidCompound> pair, double score) {
			GammaDistribution dist = new GammaDistribution(alpha, lambda);
			return dist.density(score + beta);
		}

	}

	public void find(Results results, float maxDistance, PrintWriter pw) {

		Collections.shuffle(results.getData()); // shuffle!

		for (Result result : results.getData()) {

			if (!SignificanceFactory.rotationallySymmetricSmart().isSignificant(result)) continue;

			if (result.getOrder() == null || result.getOrder() != 2) continue;

			Protodomain wholeAligned, protodomain1, protodomain2;
			
			if (result.getProtodomain() == null) continue;

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

				// k = 0.058, λ = 0.279
				SubstitutionMatrixScorer<ProteinSequence, AminoAcidCompound> scorer = new SubstitutionMatrixScorer<ProteinSequence, AminoAcidCompound>(pair, SubstitutionMatrixHelper.getBlosum62());
				double score = scorer.getSimilarity();
				logger.debug("Score: " + String.format("%.4f", score));
				
				if (score > maxDistance) continue;
				
				HmmerResult hmmer = getPfam(seq1);
				if (hmmer == null) continue;
				
				pw.println("alignment p-value: " + String.format("%.4f", score));
				pw.println("pfam name: " + hmmer.getName());
				pw.println("pfam p-value: " + String.format("%.4f", hmmer.getPvalue()));
				
				pw.println(result.getScopId());
				pw.println(protodomain1);
				pw.println(protodomain2);

				pw.println(pair.toString(80));

				pw.println("-----------------------------------------------------");

			} catch (Exception e) {
				logger.error("Encountered an error on " + result.getScopId() + " with protodomain " + result.getProtodomain(), e);
			} finally {
			}
		}
	}

}
