package org.biojava3.structure.align.symm.census;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;

import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.util.AtomCache;

public class AlignProtodomains {
	public static void main(String[] args){
		String path =  "/Volumes/Macintosh HD2/PDB/";

		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,path);


		AlignProtodomains me = new AlignProtodomains();

		me.run();

	}


	/** run protodomain searches against the set of known protodomains
	 * 
	 */
	private void run() {
		AtomCache cache = new AtomCache();

		String filePath = cache.getPath() + File.separator + "scopCensus.xml";
		List<Future<ProtoDomainDBSearchResults>> allFutureResults = new ArrayList<Future<ProtoDomainDBSearchResults>>();
		try {
			CensusResults census = ScanSCOPForSymmetry.getExistingResults(filePath);

			List<CensusResult> data = census.getData();
			
			for (CensusResult d : data){
				
				if ( ! d.getIsSignificant()) {
					continue;
				}
				
				String protoDomain = d.getProtoDomain();
				System.out.println(protoDomain);
				CallableProtodomainComparison comp = new CallableProtodomainComparison();
				comp.setCache(cache);
				comp.setData(data);
				comp.setProtoDomain(protoDomain);
				
				Future<ProtoDomainDBSearchResults> futureResults =ScanSCOPForSymmetry.pool.submit(comp);
				// run with only 1 thread...
				//futureResults.get();
				allFutureResults.add(futureResults);
			}
			
			int counter = 0;
			
			for ( Future<ProtoDomainDBSearchResults> futureResult : allFutureResults){
				
				ProtoDomainDBSearchResults pddsr = futureResult.get();
				if ( pddsr == null) {
					
					counter ++;
					continue;
				}
				
				counter ++;
				System.out.println("#" + counter + " FINISHED: " + pddsr.getProtoDomain());
			}
			
		} catch (Exception e){
			e.printStackTrace();
		}
		System.exit(0);

	}
}
