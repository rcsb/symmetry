package org.biojava3.structure.quaternary.misc;

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
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.biojava3.structure.quaternary.core.QuatSymmetryDetector;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.QuatSymmetryResults;

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

//		BlastClustReader reader100 = new BlastClustReader(100);
//		BlastClustReader reader95 = new BlastClustReader(95);
//		BlastClustReader reader30 = new BlastClustReader(30);

		//		reps.add("1AJC");
		// reps.add("1A98"); // case where seq. alignment != structural alignment due to alignment problems at gap.
		//		reps.add("2CD5"); // no bioassembly file available
		//		reps.add("4HHB"); // 2 alpha, 2 beta
		//		reps.add("2BG9"); // acetylcholin receptor, 2 alpha, 1 beta, 1 delta, 1 gamma
		//		reps.add("2WRN"); // ribosome
		//		reps.add("3SYW"); // DNA, no protein chains

		PrintWriter out = null;
//		PrintWriter out1 = null;
		PrintWriter error = null;

		try {
			out = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_symm.csv"));
//			out1 = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_error.csv"));
			error = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}


		long t1 = System.nanoTime();

		int success = 0;
		int failure = 0;

		String header = "pdbId,bioassembly,local,pseudostoichiometric,stoichiometry,pseudosymmetric,pointgroup,order,minidentity,maxidentity,rmsd,subunits, time";
		out.println(header);
//		out1.println(header);
		
		QuatSymmetryParameters parameters = new QuatSymmetryParameters();

		System.out.println("Getting PdbEntryInfo");

		//		List<PdbEntryInfo> list = PdbEntryInfoParser.getPdbEntryInfo();

		boolean skip = false;
		String restartId = "3K1P";
		Set<String> set = GetRepresentatives.getAll();

		for (String pdbId: set) {
			if (skip && pdbId.equals(restartId)) {
				skip = false;
			} 

			if (skip) {
				continue;
			}

			//			if (!pdbId.equals("1M5Q")) continue; // good example

			System.out.println("------------- " + pdbId  + "-------------");

			StructureIO.setAtomCache(cache);
			int bioAssemblyCount = StructureIO.getNrBiologicalAssemblies(pdbId);
			int bioAssemblyId = 0;
			System.out.println("Bioassemblies: " + bioAssemblyCount);
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
				} catch (StructureException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					error.println(pdbId + "[" + bioAssemblyId + "]: " + e.getMessage());
				}
				error.flush();

				long ts1 = System.nanoTime(); 	

				try {
					QuatSymmetryDetector detector = new QuatSymmetryDetector(structure, parameters);

					if (detector.hasProteinSubunits()) {	
						long ts2 = System.nanoTime();
						int time = Math.round((float)(ts2-ts1)/1000000.0f);
						for (QuatSymmetryResults results: detector.getGlobalSymmetry()) {
							out.println("PDB" + pdbId +"," + bioAssemblyId + "," + results.isLocal() +
									"," + results.getSubunits().isPseudoStoichiometric() +
									"," + results.getSubunits().getStoichiometry() +
									"," + results.getSubunits().isPseudoSymmetric() +
									"," + results.getRotationGroup().getPointGroup() +
									"," + results.getRotationGroup().getOrder() + 
									"," + Math.round(results.getSubunits().getMinSequenceIdentity()*100.0) +
									"," + Math.round(results.getSubunits().getMaxSequenceIdentity()*100.0) +
									"," + (float) results.getRotationGroup().getAverageTraceRmsd() +
									"," + results.getSubunits().getSubunitCount() +
									"," + time);
						}
						for (List<QuatSymmetryResults> resultsList: detector.getLocalSymmetries()) {
							for (QuatSymmetryResults results: resultsList) {
								out.println("PDB" + pdbId +"," + bioAssemblyId + "," + results.isLocal() +
										"," + results.getSubunits().isPseudoStoichiometric() +
										"," + results.getSubunits().getStoichiometry() +
										"," + results.getSubunits().isPseudoSymmetric() +
										"," + results.getRotationGroup().getPointGroup() +
										"," + results.getRotationGroup().getOrder() + 
										"," + Math.round(results.getSubunits().getMinSequenceIdentity()*100.0) +
										"," + Math.round(results.getSubunits().getMaxSequenceIdentity()*100.0) +
										"," + (float) results.getRotationGroup().getAverageTraceRmsd() +
										"," + results.getSubunits().getSubunitCount() +
										"," + time);
							}
						}
					}
					success++;
					out.flush();
				} catch (Exception e) {
					failure++;
					e.printStackTrace();
					error.println(pdbId + "[" + bioAssemblyId + "]: " + e.getMessage());
				}
	
				

//                ProteinComplexSignature s100 = new ProteinComplexSignature(pdbId, chainIds, reader100);
//				String signature100 = s100.getComplexSignature();
//				String stoich100 = s100.getComplexStoichiometry();
//				
//				ProteinComplexSignature s95 = new ProteinComplexSignature(pdbId, chainIds, reader95);
//				String signature95 = s95.getComplexSignature();
//				String stoich95 = s95.getComplexStoichiometry();
//				
//				ProteinComplexSignature s30 = new ProteinComplexSignature(pdbId, chainIds, reader30);
//				String signature30 = s30.getComplexSignature();
//				String stoich30 = s30.getComplexStoichiometry();
//				


		}
		long t2 = System.nanoTime();


		System.out.println("PDBs succeeded: " + success);
		System.out.println("PDBs failed   : " + failure);
		System.out.println("Total structure: " + set.size());
		System.out.println("Cpu time: " + (t2-t1)/1000000 + " ms.");
		
	    out.close();
//		out1.flush();
		error.close();
	}
	
	private void initializeCache() {
		cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(true);
		params.setLoadChemCompInfo(true);
	}
	
//	private static Structure  readStructure(String pdbId, int bioAssemblyId) {
//		// initialize the PDB_DIR env variable
//		AtomCache cache = new AtomCache();
//
//		FileParsingParameters p = new FileParsingParameters();
//		p.setStoreEmptySeqRes(true);
//		p.setLoadChemCompInfo(true);
//		p.setParseCAOnly(true);
//		p.setAtomCaThreshold(Integer.MAX_VALUE);
//		p.setParseBioAssembly(true);
//
//
//		PDBFileReader pdbreader = new PDBFileReader();
//		pdbreader.setPath(cache.getPath());
//		pdbreader.setFileParsingParameters(p);
//		pdbreader.setAutoFetch(true);
//		pdbreader.setBioAssemblyId(bioAssemblyId);
//		pdbreader.setBioAssemblyFallback(false);
//		Structure structure = null;
//		try { 
//			structure = pdbreader.getStructureById(pdbId);
//			if ( bioAssemblyId > 0 )
//				structure.setBiologicalAssembly(true);
//			structure.setPDBCode(pdbId);
//		} catch (Exception e){
//			e.printStackTrace();
//			System.exit(-1);
//		}
//		return structure;
//	}

}
