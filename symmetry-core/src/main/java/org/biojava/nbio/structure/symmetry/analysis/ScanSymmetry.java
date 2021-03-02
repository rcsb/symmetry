/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.symmetry.analysis;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.StructureFiletype;
import org.biojava.nbio.structure.chem.AllChemCompProvider;
import org.biojava.nbio.structure.chem.ChemCompGroupFactory;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetrySubunits;
import org.biojava.nbio.structure.symmetry.misc.ProteinComplexSignature;
import org.biojava.nbio.structure.symmetry.utils.BlastClustReader;
import org.biojava.nbio.structure.xtal.SpaceGroup;

/**
 * Utility to identify symmetry in a large-scale search. This is similar to
 * QuatSymmMain but less user friendly.
 *
 * usage: java ScanSymmetry pdblist resultdir
 *
 * Arguments: pdblist File containing one structure identifier per line. May be
 * - for stdin. resultdir Directory to store results
 *
 * Each run produces a DATE_symm.csv file with symmetry results and a
 * DATE_error.txt
 */
public class ScanSymmetry implements Runnable {
	private AtomCache cache = null;
	private String pdbList;
	private String resultDir;

	public ScanSymmetry(String pdbList, String resultDir) {
		this.pdbList = pdbList;
		this.resultDir = resultDir;
		initializeCache();
	}

	public static void main(String[] args) {
		if (args.length != 2) {
			System.err.print("usage: ScanSymmetry pdblist resultdir");
			System.exit(1);
		}
		String pdbList = args[0];
		String resultDir = args[1];

		if (!Files.exists(Paths.get(resultDir))) {
			System.err.printf("Directory does not exist: %s%n", resultDir);
			System.exit(1);
		}
		new ScanSymmetry(pdbList, resultDir).run();
	}

	@Override
	public void run() {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		System.out.println("Reading blastclust files");
		BlastClustReader reader95 = new BlastClustReader(95);
		BlastClustReader reader30 = new BlastClustReader(30);

		PrintWriter out = null;
		PrintWriter error = null;

		try {
			out = new PrintWriter(Files.newBufferedWriter(Paths.get(this.resultDir, timeStamp + "_symm.csv")));
			error = new PrintWriter(Files.newBufferedWriter(Paths.get(this.resultDir, timeStamp + "_error.txt")));
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}

		long t1 = System.nanoTime();

		int success = 0;
		int proteins = 0;
		int failure = 0;

		String header = "pdbId,bioassembly,local,pseudostoichiometric,stoichiometry,pseudosymmetric,pointgroup,order,"
				+ "lowSymmetry,minidentity,maxidentity,subunitrmsd,rmsd,tm,minrmsd,maxrmsd,mintm,maxtm,rmsdintra,tmintra,symdeviation,subunits,nucleiacids,cacount,time,signature95,stoich95,signature30,stoich30,spacegroup";
		out.println(header);
		out.flush();

		QuatSymmetryParameters parameters = new QuatSymmetryParameters();
		SubunitClustererParameters scParams = new SubunitClustererParameters();

		StructureIO.setAtomCache(cache);

		Collection<String> set;
		try {
			set = parsePDBList(pdbList);
		} catch (IOException e1) {
			System.err.printf("Error reading pdblist: %s%n", e1.getMessage());
			System.exit(1);
			return;
		}
		System.out.printf("Read %d identifiers%n", set.size());

		// set skip to true to restart calculation with a specified PDB ID
		boolean skip = false;
		String restartId = "10MH";

		for (String pdbId : set) {
			System.out.println(pdbId);
			if (skip && pdbId.equals(restartId)) {
				skip = false;
			}
			if (skip) {
				continue;
			}

			System.out.println("------------- " + pdbId + " -------------");

			List<Structure> structures = null;
			try {
				structures = StructureIO.getBiologicalAssemblies(pdbId);
			} catch (StructureException | IOException e) {
				e.printStackTrace();
				error.println(pdbId + ": " + e.getMessage());
				error.flush();
				continue;
			}

			int i = 0;
			for (Structure structure : structures) {

				// note: before biojava 5.0 refactoring, the structures without
				// bioassemblies would
				// use i=0 as the identifier for the default bioassembly (the
				// asymmetric unit). Now
				// identifier i=1 is used - JD 2016-05-17
				i++;

				if (!structure.getPolyChains().stream().anyMatch(Chain::isProtein)) {
					// not a protein assembly
					continue;
				}
				SpaceGroup spaceGroup = null;
				// float resolution = 0.0f;
				PDBCrystallographicInfo info = structure.getCrystallographicInfo();
				if (info != null) {
					spaceGroup = info.getSpaceGroup();
				}
				// PDBHeader pdbHeader = structure.getPDBHeader();
				// resolution = pdbHeader.getResolution();

				try {

					long ts1 = System.nanoTime();
					// save global symmetry results
					QuatSymmetryResults globalResults = QuatSymmetryDetector.calcGlobalSymmetry(structure, parameters,
							scParams);
					long ts2 = System.nanoTime();
					int time = Math.round((ts2 - ts1) / 1000000.0f);
					printToCsv(reader95, reader30, out, pdbId, i, time, globalResults, spaceGroup);

					// save local symmetry results
					ts1 = System.nanoTime();
					List<QuatSymmetryResults> localResults = QuatSymmetryDetector.calcLocalSymmetries(structure,
							parameters, scParams);
					ts2 = System.nanoTime();
					time = Math.round((ts2 - ts1) / 1000000.0f);

					for (QuatSymmetryResults localResult : localResults) {
						printToCsv(reader95, reader30, out, pdbId, i, time, localResult, spaceGroup);
					}
					proteins++;

					success++;
					out.flush();
				} catch (Exception e) {
					failure++;
					e.printStackTrace();
					error.println(pdbId + "[" + i + "]: " + e.getMessage());
					error.flush();
				}
			}
		}

		long t2 = System.nanoTime();

		System.out.println("PDBs succeeded: " + success);
		System.out.println("PDBs failed   : " + failure);
		System.out.println("Proteins      : " + proteins);
		System.out.println("Total structure: " + set.size());
		System.out.println("Cpu time: " + (t2 - t1) / 1000000 + " ms.");

		out.close();
		error.close();
	}

	private static List<String> parsePDBList(String pdbList) throws IOException {
		Stream<String> in;
		if (pdbList == null || pdbList.equals("-")) {
			in = new BufferedReader(new InputStreamReader(System.in)).lines();
		} else {
			in = Files.lines(Paths.get(pdbList));
		}
		return in.filter((line) -> line.length() > 0 && line.charAt(0) != '#').collect(Collectors.toList());
	}

	private void printToCsv(BlastClustReader reader95, BlastClustReader reader30, PrintWriter out, String pdbId,
			int bioAssemblyId, int time, QuatSymmetryResults results, SpaceGroup spaceGroup) {
		ProteinComplexSignature s95 = new ProteinComplexSignature(pdbId,
				new QuatSymmetrySubunits(results.getSubunitClusters()).getChainIds(), reader95);
		String signature95 = s95.getComplexSignature();
		String stoich95 = s95.getComplexStoichiometry();
		ProteinComplexSignature s30 = new ProteinComplexSignature(pdbId,
				new QuatSymmetrySubunits(results.getSubunitClusters()).getChainIds(), reader30);
		String signature30 = s30.getComplexSignature();
		String stoich30 = s30.getComplexStoichiometry();
		int order = 1;
		if (!results.getSymmetry().equals("H")) {
			order = results.getRotationGroup().getOrder();
		}

		out.println("PDB" + pdbId + "," + bioAssemblyId + "," + results.isLocal() + ","
				+ results.isPseudoStoichiometric() + "," + results.getStoichiometry() + "," + results.getSymmetry()
				+ "," + order + "," + isLowSymmetry(results) + "," + Math.round(100)// results.getSubunits().getMinSequenceIdentity()
																					// * 100.0)
				+ "," + Math.round(100) // results.getSubunits().getMaxSequenceIdentity()
										// * 100.0)
				+ "," + (float) results.getScores().getRmsdCenters() + "," + (float) results.getScores().getRmsd() + ","
				+ (float) results.getScores().getTm() + "," + (float) results.getScores().getMinRmsd() + ","
				+ (float) results.getScores().getMaxRmsd() + "," + (float) results.getScores().getMinTm() + ","
				+ (float) results.getScores().getMaxTm() + "," + (float) results.getScores().getRmsdIntra() + ","
				+ (float) results.getScores().getTmIntra() + "," + (float) results.getScores().getSymDeviation() + ","
				+ results.getSubunitClusters().size() + ",NA" + ","
				+ new QuatSymmetrySubunits(results.getSubunitClusters()).getCalphaCount() + "," + time + ","
				+ signature95 + "," + stoich95 + "," + signature30 + "," + stoich30 + "," + spaceGroup);
		out.flush();
	}

	private boolean isLowSymmetry(QuatSymmetryResults results) {
		return getMinFold(new QuatSymmetrySubunits(results.getSubunitClusters())) > 1
				&& results.getRotationGroup() != null && results.getRotationGroup().getPointGroup().equals("C1");
	}

	private int getMinFold(QuatSymmetrySubunits subunits) {
		if (subunits.getFolds().size() > 1) {
			return subunits.getFolds().get(1);
		}
		return subunits.getFolds().get(0);
	}

	private void initializeCache() {
		cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		cache.setFiletype(StructureFiletype.CIF);
		params.setParseCAOnly(true);
		// MmCifBiolAssemblyProvider mmcifProvider = new
		// MmCifBiolAssemblyProvider();
		// BioUnitDataProviderFactory.setBioUnitDataProvider(mmcifProvider.getClass().getCanonicalName());
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
		// ChemCompGroupFactory.setChemCompProvider(new
		// DownloadChemCompProvider());
		ChemCompGroupFactory.getChemComp("HEM"); // Force download to start
	}
}
