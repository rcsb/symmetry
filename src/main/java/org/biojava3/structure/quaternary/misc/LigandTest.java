package org.biojava3.structure.quaternary.misc;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.biojava3.structure.quaternary.core.QuatSymmetryDetector;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.QuatSymmetryResults;
import org.biojava3.structure.quaternary.core.Subunits;

public class LigandTest implements Runnable {
//	private static String PDB_PATH = "C:/Users/Peter/Documents/PDB/";
	private AtomCache cache = null;



	public LigandTest () {
		initializeCache();
	}

	public static void main(String[] args) {
		new LigandTest().run();
	}

	public void run() {

		Set<String> set = GetRepresentatives.getAll();

		// set skip to true to restart calculation with a specified PDB ID
		boolean skip = true;
		String pdbId = "1STP";

		
	
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
			} catch (StructureException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			for (int i = 0; i < structure.nrModels(); i++) {
				for (Chain chain: structure.getModel(i)) {
					for (Group g: chain.getSeqResGroups()) {
						System.out.println("group: " + g + "name: " + g.getPDBName());
				
					}
				}
			}




	
	}
	
	private void initializeCache() {
		cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
//		params.setParseCAOnly(true);
		params.setLoadChemCompInfo(true);
	}
}
