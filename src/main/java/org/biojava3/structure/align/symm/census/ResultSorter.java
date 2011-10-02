package org.biojava3.structure.align.symm.census;

import java.io.File;
import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.util.AtomCache;

public class ResultSorter {
	public static void main(String[] args){

		String path =  "/Volumes/Macintosh HD2/PDB/";

		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,path);


		AtomCache cache  = new AtomCache();

		String r = cache.getPath() + File.separator + "scopCensus.xml";

		try {
			CensusResults census = ScanSCOPForSymmetry.getExistingResults(r);

			BeanComparator zScoreComp = new BeanComparator(CensusResult.class, "getzScore", false);
			List<CensusResult> data = census.getData();
			Collections.sort(data, zScoreComp) ;
			
			System.out.println(data.get(0));
			
			BeanComparator rmsdComp = new BeanComparator(CensusResult.class, "getRmsd", true);
			
			Collections.sort(data, rmsdComp) ;
			
			System.out.println(data.get(0));
			
		
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
