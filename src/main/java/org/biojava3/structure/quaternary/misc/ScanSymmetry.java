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
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.QuatSymmetryResults;
import org.biojava3.structure.quaternary.core.QuatSymmetryDetector;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.core.Subunits;
import org.biojava3.structure.quaternary.jmolScript.JmolSymmetryScriptGenerator;

public class ScanSymmetry implements Runnable {
//	private static String PDB_PATH = "C:/Users/Peter/Documents/PDB/";
	private AtomCache cache = null;
	private static String RESULT_DIR = "C:/Users/Peter/Documents/QuatStructureComparison/";
	private static final double SEQUENCE_IDENTITY_THRESHOLD = 0.95;

	public ScanSymmetry () {
	}
	
	public static void main(String[] args) {
		new ScanSymmetry().run();
	}

	public void run() {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
	

		System.out.println("Reading blastclust files");

		BlastClustReader reader100 = new BlastClustReader(100);
		BlastClustReader reader95 = new BlastClustReader(95);
		BlastClustReader reader30 = new BlastClustReader(30);

		//		reps.add("1AJC");
		// reps.add("1A98"); // case where seq. alignment != structural alignment due to alignment problems at gap.
		//		reps.add("2CD5"); // no bioassembly file available
		//		reps.add("4HHB"); // 2 alpha, 2 beta
		//		reps.add("2BG9"); // acetylcholin receptor, 2 alpha, 1 beta, 1 delta, 1 gamma
		//		reps.add("2WRN"); // ribosome
		//		reps.add("3SYW"); // DNA, no protein chains
		
	

//	    boolean writeFile = true;
		boolean writeFile = false;
		PrintWriter out = null;
		PrintWriter out1 = null;
		PrintWriter error = null;
		
		int seqId = (int)(SEQUENCE_IDENTITY_THRESHOLD * 100);
		try {
			out = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_symm" + seqId + ".csv"));
			out1 = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_error" + seqId + ".csv"));
			error = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_error" + seqId + ".txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}


		long t1 = System.nanoTime();
		long symTime = 0;

		int multimer = 0;
		int excluded = 0;
		int total = 0;
		int err = 0;

		String header = "pdbId,bioassembly,formula,signature100,stoichiometry100,signature95,stoichiometry95,signature30,stoichiometry30,pointgroup,symops,cacount,chains,rmsdS,rmsdT,time,rx,ry,rz,jmol";
		out.println(header);
		out1.println(header);

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
			System.out.println("Bioassemblies: " + bioAssemblyCount);
			for (int i = 0; i < bioAssemblyCount; i++) {		
	
				try {
					Structure structure = StructureIO.getBiologicalAssembly(pdbId, i);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (StructureException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				long tc1 = System.nanoTime(); 	

				QuatSymmetryParameters params = new QuatSymmetryParameters();

	
//				try {
//					workflow = new QuatSymmetryDetector(structure, params);
//					if (! workflow.hasProteinSubunits()) {
//						continue;
//					}
//					QuatSymmetryResults results = workflow.getGlobalSymmetry().get(0);
//
//					RotationGroup rotationGroup = results.getRotationGroup();	
//					pointGroup = rotationGroup.getPointGroup();
//					System.out.println("Point group: " + pointGroup);
//
//					formula = results.getSubunits().getStoichiometry();
//					System.out.println("Formula: " + formula);
//
//
//					// get metrics
//					caCount = results.getSubunits().getCalphaCount();
//					groupComplete = rotationGroup.isComplete();
//					rmsd = (float) rotationGroup.getAverageSubunitRmsd();
//					rmsdT = (float) rotationGroup.getAverageTraceRmsd();
//
//					order = rotationGroup.getOrder();
//
//					Subunits subunits = results.getSubunits();
//					chainIds = results.getSubunits().getChainIds();
//					chainCount = subunits.getCenters().size();
//					AxisTransformation at = new AxisTransformation(results);
//					rx = (float) at.getDimension().x;
//					ry = (float) at.getDimension().y;
//					rz = (float) at.getDimension().z;
//					JmolSymmetryScriptGenerator g = JmolSymmetryScriptGenerator.getInstance(at,"g");
//					jmolTransform = g.getDefaultOrientation();
//				} catch (Exception e) {
//					e.printStackTrace();
//					error.println(pdbId + "------------------------------------");
//					error.println("Error during calculation: " + e.getMessage());
//					error.flush();
//					err++;
//					continue;
//				}
//				
//				
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
//				// TODO use chain signatures to label interactions of ligands
//				long tc2 = System.nanoTime();
//				long time = (tc2 - tc1)/1000000;
//				symTime += time;
//				
//				// write .csv summary file
//				if (groupComplete) {			
//					out.print(pdbId + "," + bioassemblyId + "," + formula + "," + 
//							signature100 + "," + stoich100 + "," + signature95 + "," + stoich95 + ","  + signature30 + "," + stoich30 + "," + 
//							pointGroup + "," + order + "," + caCount + "," + chainCount + "," + rmsd + "," + rmsdT + "," +
//							time + "," + rx + ","  + ry + "," + rz + "," + "\"" + jmolTransform + "\"");
//					out.println();
//					out.flush();
//					if (chainCount > 1) {
//						multimer++;
//					}
//				} else {
//					out1.print(pdbId + "," + bioassemblyId + "," + formula + "," + 
//							signature100 + "," + stoich100 + "," + signature95 + "," + stoich95 + "," + signature30 + "," + stoich30  + "," + 
//							pointGroup  + "," + order + "," + caCount + "," + chainCount  + "," + rmsd + "," + rmsdT + "," +
//							time + "," + rx + ","  + ry + "," + rz + "," +  "\"" + jmolTransform +  "\"");
//					out1.println();
//					out1.flush();
//					excluded++;
//				}
//
			}
		}
		long t2 = System.nanoTime();

		System.out.println("Cpu time: " + (t2-t1)/1000000 + " ms.");
//		System.out.println("Calc. time: " + symTime);
//		System.out.println("Multimers: " + multimer + "out of " + total);
//		System.out.println("Excluded: " + excluded);
//		System.out.println("Errors: " + err);
//		System.out.println("Total structure: " + PdbEntryInfoParser.getPdbEntryInfo().size());
//		out.close();
//		out1.flush();
//		out1.close();
//		error.close();
	}
	
	private void initializeCache() {
		cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(true);
		params.setLoadChemCompInfo(true);
	}
	
	private static Structure  readStructure(String pdbId, int bioAssemblyId) {
		// initialize the PDB_DIR env variable
		AtomCache cache = new AtomCache();

		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		p.setLoadChemCompInfo(true);
		p.setParseCAOnly(true);
		p.setAtomCaThreshold(Integer.MAX_VALUE);
		p.setParseBioAssembly(true);


		PDBFileReader pdbreader = new PDBFileReader();
		pdbreader.setPath(cache.getPath());
		pdbreader.setFileParsingParameters(p);
		pdbreader.setAutoFetch(true);
		pdbreader.setBioAssemblyId(bioAssemblyId);
		pdbreader.setBioAssemblyFallback(false);
		Structure structure = null;
		try { 
			structure = pdbreader.getStructureById(pdbId);
			if ( bioAssemblyId > 0 )
				structure.setBiologicalAssembly(true);
			structure.setPDBCode(pdbId);
		} catch (Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
		return structure;
	}

}
