package org.biojava3.structure.align.symm.census3.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.exceptions.CompoundNotFoundException;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusResultList;
import org.biojava3.structure.align.symm.census3.CensusSignificance;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;
import org.biojava3.structure.align.symm.census3.stats.CensusStatUtils;
import org.biojava3.structure.align.symm.protodomain.Protodomain;
import org.biojava3.structure.align.symm.protodomain.ProtodomainCreationException;
import org.biojava3.ws.hmmer.HmmerResult;
import org.biojava3.ws.hmmer.HmmerScan;
import org.biojava3.ws.hmmer.RemoteHmmerScan;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Find sequence motifs that correspond to symmetry subunits or protodomains.
 * @author dmyersturnbull
 */
public class ConservedSequenceFinder {

	private final static Logger logger = LoggerFactory.getLogger(ConservedSequenceFinder.class);

	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 4) {
			System.err.println("Usage: " + ConservedSequenceFinder.class.getSimpleName() + " input-census-file.xml output-file [min-score]");
			return;
		}
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(args[1]))));
		int robustnessIterations = 20;
		if (args.length > 2) {
			robustnessIterations = Integer.parseInt(args[2]);
		}
		float score = (float) 1E-4;
		if (args.length > 3) {
			score = Float.parseFloat(args[3]);
		}
		ConservedSequenceFinder finder = new ConservedSequenceFinder(new AtomCache());
		finder.find(CensusResultList.fromXML(new File(args[0])), score, pw, robustnessIterations);
		pw.close();
	}

	private AtomCache cache;

	public ConservedSequenceFinder(AtomCache cache) {
		super();
		this.cache = cache;
	}

	/**
	 * 
	 * @param alignedUnit
	 * @param order
	 * @param scopId
	 * @param score
	 * @return
	 * @throws IOException
	 * @throws ProtodomainCreationException
	 * @throws StructureException
	 * @throws CompoundNotFoundException If the Protein structure contained an unrecognized amino acid
	 */
	public List<String> getConservedPfamNames(String alignedUnit, int order, String scopId, float score) throws IOException, ProtodomainCreationException, StructureException, CompoundNotFoundException {

		Protodomain[]  protodomains;
		Atom[][] ca;
		List<ProteinSequence> seqs;

		Protodomain wholeAligned = Protodomain.fromString(alignedUnit, scopId, cache);
		protodomains = new Protodomain[order];
		ca = new Atom[order][];
		seqs = new ArrayList<ProteinSequence>(order);
		for (int i = 0; i < order; i++) {
			protodomains[i] = wholeAligned.createSubstruct(order, i);
			ca[i] = cache.getAtoms(protodomains[i].getString());
			ProteinSequence seq = new ProteinSequence(StructureTools.convertAtomsToSeq(ca[i]));
			seqs.add(seq);
		}

		List<String> pfamNames = getConservedPfamNames(seqs, score);
		return pfamNames;
	}

	public List<String> getConservedPfamNames(List<ProteinSequence> seqs, double eValue) throws IOException {
		HmmerScan hmmer = new RemoteHmmerScan();
		Map<String, Integer> counts = new HashMap<String, Integer>();
		for (ProteinSequence seq : seqs) {
			SortedSet<HmmerResult> results = hmmer.scan(seq);
			logger.debug("Found " + results.size() + " HMMER results");
			for (HmmerResult result : results) {
				logger.debug("Found result " + result);
				if (result.getEvalue() <= eValue) {
					CensusStatUtils.plus(counts, result.getName());
				}
			}
		}
		List<String> conserved = new ArrayList<String>(1);
		for (Map.Entry<String, Integer> entry : counts.entrySet()) {
			if (entry.getValue() == seqs.size()) conserved.add(entry.getKey());
		}
		return conserved;
	}

	public static CensusResultList get1PerSuperfamily(CensusResultList results, int robustnessIterations) {
		Set<String> scopIds = new HashSet<String>();
		CensusResultList filtered = new CensusResultList();
		for (int iter = 0; iter < robustnessIterations; iter++) {
			filtered.setMeanSecondsTaken(results.getMeanSecondsTaken());
			filtered.setStartingTime(results.getStartingTime());
			Collections.shuffle(results.getEntries()); // shuffle!
			ScopDatabase scop = ScopFactory.getSCOP();
			Set<Integer> sfs = new HashSet<Integer>();
			for (CensusResult result : results.getEntries()) {
				if (scopIds.contains(result.getId())) continue;
				scopIds.add(result.getId());
				int sf = scop.getDomainByScopID(result.getId()).getSuperfamilyId();
				if (!sfs.contains(sf)) {
					filtered.add(result);
					sfs.add(sf);
				}
			}
		}
		return filtered;
	}

	public void find(CensusResultList originalResults, float score, PrintWriter pw, int robustnessIterations) {

		CensusResultList results = get1PerSuperfamily(originalResults, robustnessIterations);

		CensusSignificance sig = CensusSignificanceFactory.forCeSymmOrd();

		for (CensusResult result : results.getEntries()) {

			if (!sig.isSignificant(result)) continue;

			if (result.getGroup() == null || result.getGroup().isAsymmetric()) continue;

			if (result.getAlignedUnit() == null) continue;

			try {

				List<String> pfamNames = getConservedPfamNames(result.getAlignedUnit(), result.getOrder(), result.getId(), score);

				if (pfamNames.isEmpty()) continue;
				logger.debug("Found " + pfamNames.size() + " pfam hits for " + result.getAlignedUnit());

				pw.print(result.getId() + "\t" + result.getAlignedUnit());
				for (String name : pfamNames) {
					pw.print("\t" + name);
				}
				pw.println();
				pw.flush();

			} catch (Exception e) {
				logger.error("Encountered an error on " + result.getId() + " with protodomain " + result.getAlignedUnit(), e);
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
