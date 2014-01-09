package org.biojava3.structure.align.symm.census2.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;
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
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(args[1]))));
		double score = Math.pow(10, -10);
		if (args.length > 2) {
			score = Float.parseFloat(args[1]);
		}
		ConservedSequenceFinder finder = new ConservedSequenceFinder(new AtomCache());
		finder.find(Results.fromXML(new File(args[0])), score, pw);
		pw.close();
	}

	private AtomCache cache;

	public ConservedSequenceFinder(AtomCache cache) {
		super();
		this.cache = cache;
	}

	private String getConservedPfamName(List<ProteinSequence> seqs, double eValue) throws IOException {
		HmmerScan hmmer = new RemoteHmmerScan();
		Map<String, Integer> counts = new HashMap<String, Integer>();
		for (ProteinSequence seq : seqs) {
			SortedSet<HmmerResult> results = hmmer.scan(seq);
			logger.debug("Found " + results.size() + " HMMER results");
			for (HmmerResult result : results) {
				logger.debug("Found result " + result);
				if (result.getEvalue() <= eValue) {
					StatUtils.plus(counts, result.getName());
				}
			}
		}
		for (Map.Entry<String, Integer> entry : counts.entrySet()) {
			if (entry.getValue() == seqs.size()) return entry.getKey();
		}
		return null;
	}

	public void find(Results results, double score, PrintWriter pw) {

		Collections.shuffle(results.getData()); // shuffle!

		for (Result result : results.getData()) {

			if (!SignificanceFactory.rotationallySymmetricSmart().isSignificant(result)) continue;

			if (result.getOrder() == null || result.getOrder() < 2) continue;

			int order = result.getOrder();

			Protodomain wholeAligned;

			Protodomain[]  protodomains;
			Atom[][] ca;
			List<ProteinSequence> seqs;

			if (result.getProtodomain() == null) continue;

			try {

				wholeAligned = Protodomain.fromString(result.getProtodomain(), result.getScopId(), cache);
				protodomains = new Protodomain[order];
				ca = new Atom[order][];
				seqs = new ArrayList<ProteinSequence>(order);
				for (int i = 0; i < order; i++) {
					protodomains[i] = wholeAligned.createSubstruct(result.getOrder(), i).spliceApproxConsecutive();
					protodomains[i].buildStructure();
					ca[i] = StructureTools.getAtomCAArray(protodomains[i].getStructure());
					seqs.set(i, new ProteinSequence(StructureTools.convertAtomsToSeq(ca[i])));
				}

				//				Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(seqs);
				//				double score = scoreProfile(profile, SubstitutionMatrixHelper.getBlosum62());
				//				logger.debug("Score: " + String.format("%.4f", score));
				//				if (score > maxDistance) continue;

				String pfamName = getConservedPfamName(seqs, score);
				if (pfamName == null) continue;

				//				pw.println("alignment score: " + String.format("%5.8e", score));
				pw.println(result.getScopId() + ": \t" + pfamName);

				pw.println(result.getScopId());
				for (Protodomain protodomain : protodomains) {
					pw.println(protodomain);
				}

//				pw.println(profile.toString(80));
				pw.println("-----------------------------------------------------");

			} catch (Exception e) {
				logger.error("Encountered an error on " + result.getScopId() + " with protodomain " + result.getProtodomain(), e);
			} finally {
			}
		}
	}

	private <S extends Sequence<C>, C extends Compound> double scoreProfile(Profile<S,C> profile, SubstitutionMatrix<C> substitutionMatrix) {

		// make a RealMatrix copy of substitutionMatrix
		RealMatrix matrix = new Array2DRowRealMatrix(substitutionMatrix.getMatrix().length, substitutionMatrix.getMatrix()[0].length);
		for (int row = 0; row < matrix.getRowDimension(); row++) {
			for (int column = 0; column < matrix.getColumnDimension(); column++) {
				matrix.setEntry(row, column, substitutionMatrix.getMatrix()[row][column]);
			}
		}

		// calculate the sum of the 2-norms
		RealVector sumVector = new ArrayRealVector(profile.getLength());
		for (int column = 1; column <= profile.getLength(); column++) {
			
			// make a vector of the counts
			int[] intCounts = profile.getCompoundCountsAt(column);
			double[] counts = new double[intCounts.length];
			for (int j = 0; j < intCounts.length; j++) {
				counts[j] = intCounts[j];
			}
			
			// now record ||matrix * vector||
			RealVector vector = new ArrayRealVector(matrix.operate(counts));
			sumVector.setEntry(column, vector.getNorm());
		}
		return sumVector.getLInfNorm();
		
	}

}
