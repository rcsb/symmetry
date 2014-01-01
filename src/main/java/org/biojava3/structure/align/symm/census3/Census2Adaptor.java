package org.biojava3.structure.align.symm.census3;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava3.structure.align.symm.census2.Alignment;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;

/**
 * Converts a census2 XML file to a census3 one.
 * @author dmyersturnbull
 */
public class Census2Adaptor {

	private static final Logger logger = LogManager.getLogger(Census2Adaptor.class.getSimpleName());

	public static void main(String[] args) throws IOException {
		
		if (args.length != 2) {
			System.err.println("Usage: " + Census2Adaptor.class.getSimpleName() + " input-census2-file.xml output-census3-dir");
			return;
		}
		
		logger.info("Reading in...");
		Results census2 = Results.fromXML(args[0]);
		AtomCache cache = new AtomCache();
		
		convertInChunks(census2, cache, args[1], 1000);
//		CensusResultList newResults = convertResults(census2, cache);
//
//		// print the new census
//		IOUtils.print(newResults.toXML(), new File(args[1]));
	}

	public static void convertInChunks(Results census2, AtomCache cache, String outputDir, int every) throws IOException {

		if (!outputDir.endsWith(File.separator)) outputDir += File.separator;
		
		// conserve memory
		logger.info("Discarding unnecessary fields...");
		for (Result result : census2.getData()) {
			result.setDescription(null);
			result.setClassification(null);
			result.setSunId(null);
			result.setIsSignificant(null);
			result.setRank(null);
			result.getAlignment().setAlignScore(null);
			result.getAlignment().setCoverage(null);
			result.getAlignment().setBlock1Length(null);
			result.getAlignment().setBlock2Length(null);
		}

		CensusResultList newResults = makeEmptyResultList(census2);

		int start = 0;
//		while (getFile(outputDir, every, start, 1).exists()) {
//			start += every;
//		}
		File file = getFile(outputDir, every, start, 0);
		
		for (int i = start; i < census2.size(); i++) {
			
			Result result = census2.getData().get(i);

			try {
				CensusResult newResult = convertResult(result, cache);
				newResults.add(newResult);
			} catch (RuntimeException e) {
				logger.error("Encountered an error for " + result.getScopId(), e);
			}
			
			if (i > start && i % every == 0) {
				logger.info("Printing " + every + " results to " + file.getPath());
				print(newResults.toXML(), file);
				newResults = makeEmptyResultList(census2);
//				while (getFile(outputDir, every, i, 1).exists()) {
//					i += every;
//				}
				file = getFile(outputDir, every, i, 0);
			}
		}
		
		file = new File(outputDir + ((census2.size() / every) + 1) + ".xml");
		print(newResults.toXML(), file);
	}
		
	public static CensusResultList convertResults(Results census2, AtomCache cache) {
		
		CensusResultList newResults = makeEmptyResultList(census2);
		
		for (Result result : census2.getData()) {
			try {
				CensusResult newResult = convertResult(result, cache);
				newResults.add(newResult);
			} catch (RuntimeException e) {
				logger.error("Encountered an error for " + result.getScopId(), e);
			}
		}
		
		return newResults;
	}
	
	private static CensusResultList makeEmptyResultList(Results census2) {

		CensusResultList newResults = new CensusResultList();
		
		newResults.setMeanSecondsTaken((float) (census2.getMeanSecondsTaken() / 1000.0));
		newResults.setStartingTime(census2.getTimestamp());
		
		return newResults;
	}
	
	private static CensusResult convertResult(Result result, AtomCache cache) {
		CensusResult newResult = new CensusResult();
		newResult.setId(result.getScopId());
		newResult.setAlignedUnit(result.getProtodomain());
		if (result.getOrder() != null) {
			newResult.setGroup(new CensusSymmetryGroup("C" + result.getOrder()));
		}
		newResult.setScoreList(convertScores(result.getAlignment()));
		if (result.getAlignmentMapping() != null) {
			newResult.setAlignment(new CensusAlignment(result.getAlignmentMapping().getSimpleFunction()));
			newResult.setAxis(convertAxis(newResult.getId(), newResult.getAlignment(), cache));
		}
		return newResult;
	}
	
	private static File getFile(String outputDir, int every, int i, int offset) {
		return  new File(outputDir + ((i / every) + offset) + ".xml");
	}

	private static CensusScoreList convertScores(Alignment alignment) {
		CensusScoreList list = new CensusScoreList();
		list.setAlignLength(alignment.getAlignLength());
		list.setGapLength(alignment.getGapLength());
		list.setIdentity(alignment.getIdentity());
		list.setRmsd(alignment.getRmsd());
		list.setSimilarity(alignment.getSimilarity());
		list.setTmScore(alignment.getTmScore());
		list.setzScore(alignment.getzScore());
		return list;
	}

	private static CensusAxis convertAxis(String scopId, CensusAlignment alignment, AtomCache cache) {
		try {
			Atom[] ca1 = cache.getAtoms(scopId);
			Atom[] ca2 = cache.getAtoms(scopId);
			AFPChain afpChain = alignment.buildAfpChain(ca1, ca2);
			CensusAxis axis = new CensusAxis(new RotationAxis(afpChain));
			return axis;
		} catch (Exception e) {
			logger.error("Couldn't create axis for " + scopId, e);
			return null;
		}
	}

	/**
	 * Prints {@code string} to {@code file}.
	 * @throws IOException If any I/O exception occurs while printing; this method does not catch any exceptions
	 * @deprecated Replace with IOUtils.print.
	 */
	public static void print(String string, File file) throws IOException {
		PrintWriter out = null;
		try {
			out = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			out.print(string);
			out.flush();
			out.close();
		} finally {
			if (out != null) out.close();
		}
	}

}
