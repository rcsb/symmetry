package org.biojava3.structure.align.symm.quarternary;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.dbscan.GetRepresentatives;

public class BAFileReaderTest {

	private static String PDB_PATH = "C:/PDB/";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		AtomCache cache = new AtomCache();
		cache.setAutoFetch(true);
		cache.setPath(PDB_PATH);	
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		
		p.setAtomCaThreshold(Integer.MAX_VALUE);
	//	p.setAcceptedAtomNames(new String[]{" CA ", " CB "});
		cache.setFileParsingParams(p);

		PrintWriter error = null;
		try {

			error = new PrintWriter(new FileWriter(PDB_PATH + "error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		int total = 0;
		int exceptions = 0;	
	        
		Set<String> reps = GetRepresentatives.getAll();
		
		for (String pdbId : reps){
//			
//			if (!pdbId.equals("2DAZ")) {
//				continue;
//			}
			System.out.println("------------- " + pdbId  + "-------------");
			if (pdbId.equals("1M4X")) continue; // largest PDB assembly, causes occasional GC error

			try {
				Structure structure = cache.getBiologicalAssembly(pdbId, 1, true);
				//				structure = cache. getStructure(pdbId);

			} catch (StructureException e) {
				exceptions++;
				error.println(pdbId + "------------------------------------");
				error.println(e.getStackTrace());
		//		error.println(e.getMessage());
				error.flush();
				continue;
			} catch (IOException e) {
				exceptions++;
				error.println(pdbId + "------------------------------------");
				error.println(e.getStackTrace());
				error.flush();
				continue;
			}
			
			total++;
		}
		
		error.close();
	}
}
