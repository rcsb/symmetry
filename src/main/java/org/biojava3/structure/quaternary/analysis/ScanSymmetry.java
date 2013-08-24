package org.biojava3.structure.quaternary.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.biojava3.structure.quaternary.core.AxisAligner;
import org.biojava3.structure.quaternary.core.QuatSymmetryDetector;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.QuatSymmetryResults;
import org.biojava3.structure.quaternary.core.Subunits;
import org.biojava3.structure.quaternary.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava3.structure.quaternary.misc.ProteinComplexSignature;
import org.biojava3.structure.quaternary.utils.BlastClustReader;

public class ScanSymmetry implements Runnable {
	//	private static String PDB_PATH = "C:/Users/Peter/Documents/PDB/";
	private AtomCache cache = null;
	private static String RESULT_DIR = "C:/Users/Peter/Documents/QuatStructureComparison/";


	public ScanSymmetry () {
		initializeCache();
	}

	public static void main(String[] args) {
		new ScanSymmetry().run();
	}

	public void run() {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		System.out.println("Reading blastclust files");

		BlastClustReader reader95 = new BlastClustReader(95);
		BlastClustReader reader30 = new BlastClustReader(30);


		PrintWriter out = null;
		PrintWriter error = null;

		try {
			out = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_symm.csv"));
			error = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}


		long t1 = System.nanoTime();

		int success = 0;
		int proteins = 0;
		int failure = 0;

		String header = "pdbId,bioassembly,local,pseudostoichiometric,stoichiometry,pseudosymmetric,pointgroup,order," +
				"lowSymmetry,minidentity,maxidentity,rmsd,tm,minrmsd,maxrmsd,mintm,maxtm,subunits,nucleiacids,time,signature95,stoich95,signature30,stoich30,spacegroup";
		out.println(header);

		QuatSymmetryParameters parameters = new QuatSymmetryParameters();

		parameters.setVerbose(true);

		Set<String> set = GetRepresentatives.getAll();

		// set skip to true to restart calculation with a specified PDB ID
		boolean skip = false;
		String restartId = "10MH";

		for (String pdbId: set) {
//		for (String pdbId: testCase) {
//					for (String pdbId: helixExamples) {
//		for (String pdbId: collagenExamples) {
			if (skip && pdbId.equals(restartId)) {
				skip = false;
			} 
			if (skip) {
				continue;
			}

			// exclude the following examples (out of memory exception)		
			if (pdbId.equals("1M4X")) continue;
			if (pdbId.equals("3HQV")) continue;
			if (pdbId.equals("3HR2")) continue;
			if (pdbId.equals("4A8B")) continue; 
			if (pdbId.equals("4D8Q")) continue;

			if (pdbId.equals("4A0W")) continue;

			System.out.println("------------- " + pdbId  + "-------------");

			StructureIO.setAtomCache(cache);
			int bioAssemblyCount = StructureIO.getNrBiologicalAssemblies(pdbId);

			ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

			int bioAssemblyId = 0;
			System.out.println("Bioassemblies: " + bioAssemblyCount);
			System.out.println(ChemCompGroupFactory.getChemCompProvider().getClass().getName());
			if (bioAssemblyCount > 0) {
				bioAssemblyId = 1;
			}

			System.out.println("bioAssemblyId: " + bioAssemblyId);
			//			for (int i = 0; i < bioAssemblyCount; i++) {	
			Structure structure = null;
			try {
				structure = StructureIO.getBiologicalAssembly(pdbId, bioAssemblyId);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				error.println(pdbId + "[" + bioAssemblyId + "]: " + e.getMessage());
				error.flush();
			} catch (StructureException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				error.println(pdbId + "[" + bioAssemblyId + "]: " + e.getMessage());
				error.flush();
			}

			long ts1 = System.nanoTime(); 	

			try {
				String spaceGroup = "";
				if (structure != null) {
					spaceGroup = structure.getCrystallographicInfo().getSpaceGroup();
				}
				QuatSymmetryDetector detector = new QuatSymmetryDetector(structure, parameters);

				if (detector.hasProteinSubunits()) {	
					long ts2 = System.nanoTime();

					int time = Math.round((float)(ts2-ts1)/1000000.0f);
					
					// save global symmetry results
					List<QuatSymmetryResults> globalResults = detector.getGlobalSymmetry();
					printToCsv(reader95, reader30, out, pdbId,
							bioAssemblyId, time, globalResults, spaceGroup);
					
					// save local symmetry results
					for (List<QuatSymmetryResults> localResults: detector.getLocalSymmetries()) {
						printToCsv(reader95, reader30, out, pdbId,
								bioAssemblyId, time, localResults, spaceGroup);
					}
					proteins++;
				}
				success++;
				out.flush();
			} catch (Exception e) {
				failure++;
				e.printStackTrace();
				error.println(pdbId + "[" + bioAssemblyId + "]: " + e.getMessage());
				error.flush();
			}
		}
		long t2 = System.nanoTime();

		System.out.println("PDBs succeeded: " + success);
		System.out.println("PDBs failed   : " + failure);
		System.out.println("Proteins      : " + proteins);
		System.out.println("Total structure: " + set.size());
		System.out.println("Cpu time: " + (t2-t1)/1000000 + " ms.");

		out.close();
		//		out1.flush();
		error.close();
	}

	private void printToCsv(BlastClustReader reader95,
			BlastClustReader reader30, PrintWriter out, String pdbId,
			int bioAssemblyId, int time, List<QuatSymmetryResults> resultsList, String spaceGroup) {
		for (QuatSymmetryResults results: resultsList) {
			ProteinComplexSignature s95 = new ProteinComplexSignature(pdbId, results.getSubunits().getChainIds(), reader95);
			String signature95 = s95.getComplexSignature();
			String stoich95 = s95.getComplexStoichiometry();
			ProteinComplexSignature s30 = new ProteinComplexSignature(pdbId, results.getSubunits().getChainIds(), reader30);
			String signature30 = s30.getComplexSignature();
			String stoich30 = s30.getComplexStoichiometry();
			int order = 1;
			if (!results.getSymmetry().equals("H")) {
				order = results.getRotationGroup().getOrder();
			}
			AxisAligner aligner = AxisAligner.getInstance(results);
			JmolSymmetryScriptGenerator script = JmolSymmetryScriptGenerator.getInstance(aligner, "g");
			String color = script.colorBySymmetry();
			String orient = script.getOrientationWithZoom(0);

			out.println("PDB" + pdbId +"," + bioAssemblyId + "," + results.isLocal() +
					"," + results.getSubunits().isPseudoStoichiometric() +
					"," + results.getSubunits().getStoichiometry() +
					"," + results.getSubunits().isPseudoSymmetric() +
					"," + results.getSymmetry() +
					"," + order + 
					"," + isLowSymmetry(results) +
					"," + Math.round(results.getSubunits().getMinSequenceIdentity()*100.0) +
					"," + Math.round(results.getSubunits().getMaxSequenceIdentity()*100.0) +
					"," + (float) results.getScores().getRmsd() +
					"," + (float) results.getScores().getTm() +
					"," + (float) results.getScores().getMinRmsd() +
					"," + (float) results.getScores().getMaxRmsd() +
					"," + (float) results.getScores().getMinTm() +
					"," + (float) results.getScores().getMaxTm() +
					"," + results.getSubunits().getSubunitCount() +
					"," + results.getNucleicAcidChainCount() +
					"," + time +
					"," + signature95 +
					"," + stoich95 +
					"," + signature30 +
					"," + stoich30 +
					"," + spaceGroup +
					"," + "\"" + orient +  "\"" +
					"," +  "\"" + color +  "\""
					);
		}
	}

	private boolean isLowSymmetry(QuatSymmetryResults results) {
		return getMinFold(results.getSubunits()) > 1 && results.getRotationGroup() != null && results.getRotationGroup().getPointGroup().equals("C1");
	}

	private int getMinFold(Subunits subunits) {
		if (subunits.getFolds().size() > 1) {
			return subunits.getFolds().get(1);
		}
		return subunits.getFolds().get(0);
	}

	private void initializeCache() {
		cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(true);
		params.setLoadChemCompInfo(true);
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
	}
	
	private static String[] testCase = {"1A02"};
	
	// global helix examples from ccquickly
	private static String[] helixExamples = {
	"1B47",
	"1BKV",
	"1C09",
	"1CGD",
	"1CGM",
	"1CR0",
	"1CR1",
	"1CR2",
	"1CR4",
	"1FC3",
	"1FFX",
	"1G9F",
	"1HGV",
	"1HGZ",
	"1HH0",
	"1IFD",
	"1IFI",
	"1IFJ",
	"1IFK",
	"1IFL",
	"1IFM",
	"1IFN",
	"1IFP",
	"1JI7",
	"1K6F",
	"1L5A",
	"1MM9",
	"1MOY",
	"1N03",
	"1NMT",
	"1PFI",
	"1PV4",
	"1PVO",
	"1Q7D",
	"1QL1",
	"1QL2",
	"1QSU",
	"1R6Z",
	"1RIR",
	"1RMV",
	"1RQ0",
	"1SZP",
	"1T5E",
	"1U94",
	"1U98",
	"1U99",
	"1VF7",
	"1VTM",
	"1VZJ",
	"1WUD",
	"1WZB",
	"1XMS",
	"1XMV",
	"1XP8",
	"1YJ7",
	"1YM8",
	"1YS3",
	"1YSR",
	"1Z0B",
	"1Z0C",
	"1Z0W",
	"1Z4V",
	"1Z4W",
	"1Z4X",
	"1Z4Y",
	"1Z4Z",
	"1Z50",
	"1ZKK",
	"2AVP",
	"2C0W",
	"2CUO",
	"2D3F",
	"2D3H",
	"2DRT",
	"2FKJ",
	"2FUF",
	"2G66",
	"2HCB",
	"2HI5",
	"2HIL",
	"2HY6",
	"2IEF",
	"2IFM",
	"2IFN",
	"2IFO",
	"2KJ3",
	"2LBU",
	"2LPZ",
	"2OM3",
	"2QU4",
	"2R19",
	"2R1A",
	"2TMV",
	"2V4D",
	"2V6L",
	"2VE9",
	"2WYY",
	"2XEA",
	"2Y83",
	"3A08",
	"3A0A",
	"3A0M",
	"3A19",
	"3A1H",
	"3ABN",
	"3ADM",
	"3AH9",
	"3AI6",
	"3B0S",
	"3B2C",
	"3BQ7",
	"3BYH",
	"3DMW",
	"3DTP",
	"3EDL",
	"3F4Z",
	"3G37",
	"3HP3",
	"3IFM",
	"3IKU",
	"3IKY",
	"3IPN",
	"3J06",
	"3J0R",
	"3J0S",
	"3J1R",
	"3J2U",
	"3JQH",
	"3KQU",
	"3KZE",
	"3LG7",
	"3LUE",
	"3LVH",
	"3M6A",
	"3MFP",
	"3MOP",
	"3NTU",
	"3OQ9",
	"3P46",
	"3PDM",
	"3POD",
	"3PON",
	"3PXI",
	"3QIL",
	"3R8F",
	"3RD4",
	"3T98",
	"3U29",
	"4A6J",
	"4A7N",
	"4AXY",
	"4BS1",
	"4BT0",
	"4DG7",
	"4DMT",
	"4E2H",
	"4GHA",
	"4GHL",
	"4GYX",
	"4IFM"
	};


	private static String[] helixExamples2 = {
		"1B47","1BKV","1C09","1CGD","1CGM","1CR0","1CR1","1CR2","1CR4","1FC3","1FFX","1FZD","1GL2","1HGV","1HGZ","1HH0","1IFD","1IFI",
		"1IFJ","1IFK","1IFL","1IFM","1IFN","1IFP","1JI7","1K6F","1L5A","1M8Q","1MM9","1MOY","1MVW","1MVW","1N03","1NMT","1O18",
		"1O19","1O1A","1O1B","1O1C","1O1D","1O1E","1O1F","1O1G","1PFI","1PV4","1PVO","1QL1","1QL2","1QSU","1QVR","1R6Z","1RHG","1RIR",
		"1RMV","1RQ0","1SA0","1SA1","1SZP","1T5E","1U94","1U98","1U99","1VF7","1VTM","1VZJ","1WUD","1XMS","1XMV","1XP8","1YJ7","1YS3",
		"1YSR","1Z0B","1Z0C","1Z0W","1Z2B","1Z4V","1Z4W","1Z4X","1Z4Y"};

	private static String[] collagenExamples = {
		"1A3I","1A3J","1BKV","1CAG","1CGD","1DZI","1EAK","1EI8","1G9W","1ITT",
		"1K6F","1NAY","1Q7D","1QSU","1V4F","1V6Q","1V7H","1WZB","1X1K","1YM8",
		"2CUO","2D3F","2D3H","2DRT","2DRX","2F6A","2G66","2KLW","2LLP","2V53",
		"2WUH","2Y5T","3A08","3A0A","3A0M","3A19","3A1H","3ABN","3ADM","3AH9",
		"3AI6","3B0S","3B2C","3DMW","3IPN","3P46","3POB","3POD","3PON","3T4F",
		"3U29","3ZHA","4AU2","4AU3","4AUO","4AXY","4DMT","4DMU","4GYX"
	};
}
