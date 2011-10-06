package demo;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;

import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.core.util.ConcurrencyTools;
import org.biojava3.structure.align.symm.census.CallableProtodomainComparison;
import org.biojava3.structure.align.symm.census.CensusResult;
import org.biojava3.structure.align.symm.census.CensusResults;
import org.biojava3.structure.align.symm.census.ProtoDomainDBSearchResults;
import org.biojava3.structure.align.symm.census.ScanSCOPForSymmetry;

public class AlignProtodomains {
	
	static {
		int maxThreads = Runtime.getRuntime().availableProcessors()-1;
		if ( maxThreads < 1)
			maxThreads = 1;

		System.out.println("using " + maxThreads + " threads for alignment calculation of DB jobs...");
		//pool = Executors.newFixedThreadPool(maxThreads);
		ConcurrencyTools.setThreadPoolSize(maxThreads);
	}

	
	
	public static void main(String[] args){
		String path =  "/Users/ap3/WORK/PDB/";

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
				
				Future<ProtoDomainDBSearchResults> futureResults =ConcurrencyTools.submit(comp);
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
