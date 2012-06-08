package org.biojava3.structure.align.quaternary;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.align.symm.quaternary.analysis.PdbEntryInfo;
import org.biojava3.structure.align.symm.quaternary.analysis.PdbEntryInfoParser;
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
	        
		for (PdbEntryInfo entry: PdbEntryInfoParser.getPdbEntryInfo()){
			String pdbId = entry.getPdbId();
			System.out.println("------------- " + pdbId  + "-------------");

			int bioAssemblyCount = entry.getBioAssemblyCount();
			for (int i = 0; i < bioAssemblyCount; i++) {
				if (pdbId.equals("1M4X")) continue; // largest PDB assembly, causes occasional GC error

				Structure structure = null;
				try {
					if (bioAssemblyCount == 0) {
						structure = cache. getStructure(pdbId);
					} else {
						structure = cache.getBiologicalAssembly(pdbId, i+1, true);
					}

				} catch (StructureException e) {
					e.printStackTrace();
					error.println(pdbId + "------------------------------------");
					error.println(e.getMessage());
					error.flush();
					exceptions++;
					continue;
				} catch (IOException e) {
					e.printStackTrace();
					error.println(pdbId + "------------------------------------");
					error.println(e.getMessage());
					error.flush();
					exceptions++;
					continue;
				}

			total++;
		}
	}
		
		error.println("total stuctures: " + total + " exceptions: " + exceptions);
		error.close();
	}
}
