package org.biojava3.structure.align.symm.census;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.concurrent.Future;

import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.SynchronizedOutFile;
import org.biojava.bio.structure.domain.DomainProvider;
import org.biojava.bio.structure.domain.DomainProviderFactory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.core.util.ConcurrencyTools;
import org.biojava3.structure.utils.FileUtils;
import org.rcsb.fatcat.server.PdbChainKey;

public class ScanRepresentativesForSymmetry {

	public static String newline = System.getProperty("line.separator");


	static {
		int maxThreads = Runtime.getRuntime().availableProcessors()-1;
		if ( maxThreads < 1)
			maxThreads = 1;

		System.out.println("using " + maxThreads + " threads for alignment calculation of DB jobs...");
		//pool = Executors.newFixedThreadPool(maxThreads);
		ConcurrencyTools.setThreadPoolSize(maxThreads);
	}


	public static void main(String[] args){

		String path =  "/Users/andreas/WORK/PDB/";

		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,path);

		try {

			ScanRepresentativesForSymmetry me = new ScanRepresentativesForSymmetry();

			me.run();

		} catch (Exception e){
			e.printStackTrace();
		}

	}
	private void run() {
		AtomCache cache = new AtomCache();

		DomainProvider domainProvider = DomainProviderFactory.getDomainProvider();

		SortedSet<String> representatives = domainProvider.getRepresentativeDomains();

		String f = cache.getPath() + File.separator + "scopCensus.html";

		String resultsXmlFileLocation = cache.getPath() + File.separator + "scopCensus.xml";

		System.out.println("results will be at " + resultsXmlFileLocation);
		try {

			CensusResults census = getExistingResults(resultsXmlFileLocation) ;
			if ( census == null)
				census = new CensusResults();

			List<String> knownResults = getKnownResults(census);
			int count = 0;



			List<Future<CensusResult>> futureData = new ArrayList<Future<CensusResult>>();

			ScopDatabase scop = ScopFactory.getSCOP();


			for (String name: representatives){

				PdbChainKey key = PdbChainKey.fromName(name);


				count++;
				//ScopDomain first = null;
				ScopDescription superfamily = null;
				if ( key.isScopName()) {
					ScopDomain scopDomain = scop.getDomainByScopID(name);
					int sunid = scopDomain.getSuperfamilyId();
					superfamily = scop.getScopDescriptionBySunid(sunid);
					//int sunid = superfamily.getSunID();
					//List<ScopDomain> familyMembers = scop.getScopDomainsBySunid(sunid);
					//first = familyMembers.get(0);
				}

				if ( knownResults.contains(name))
					continue;

				if ( name.equals("ds046__"))
					continue;

				SymmetryCalculationJob calc = new SymmetryCalculationJob();
				calc.setName(name);
				calc.setCache(cache);
				calc.setScopDescription(superfamily);
				calc.setCount(count);

				Future<CensusResult> result = ConcurrencyTools.submit(calc);
				futureData.add(result);
			}

			File o = new File(f);
			if ( o.exists())
				o.delete();
			SynchronizedOutFile outFile = new SynchronizedOutFile(o);


			ResultConverter.printHTMLHeader(outFile);

			//int fragmentLength = 8;

			int withSymm = 0;
			Map<Character,Integer> classStats = new HashMap<Character, Integer>();
			Map<Character,Integer> totalStats = new HashMap<Character, Integer>();

			List<CensusResult> allResults = census.getData();

			withSymm = countNrSymm(outFile, withSymm, classStats, totalStats, allResults);

			System.out.println("starting to process multi threaded results...");
			for (Future<CensusResult> futureResult : futureData) {

				CensusResult result = futureResult.get();
				System.out.println("got result: " + result);
				if ( result == null)
					continue;

				boolean isSymmetric = processResult(result, outFile,totalStats,classStats);
				if ( isSymmetric)
					withSymm++;

				allResults.add(result);
				//if ( withSymm > 5)
				//	break;
				if ( allResults.size() % 100 == 0)
					writeResults(census, resultsXmlFileLocation);
			} 

			outFile.write("</tbody></table>");

			census.setData(allResults);

			writeResults(census, resultsXmlFileLocation);


			System.out.println("===");
			System.out.println("Overall symmetry: " + String.format("%.2f",(withSymm/(float)count)) + "%" );
			System.out.println("Statistics for SCOP classes:");
			for (Character scopClass: totalStats.keySet()){
				Integer total = totalStats.get(scopClass);
				Integer symm  = classStats.get(scopClass);
				System.out.println(scopClass);
				System.out.println("Class: " + scopClass + " " + String.format("%.2f",(symm/(float)total))  + "%");
			}

			outFile.close();
		} catch (Exception e){
			e.printStackTrace();
		}




	}


	private int countNrSymm(SynchronizedOutFile outFile, int withSymm,
			Map<Character, Integer> classStats,
			Map<Character, Integer> totalStats, List<CensusResult> allResults)
					throws IOException {
		// process know results again to make sure the stats are correct
		for ( CensusResult result : allResults){
			boolean isSymmetric = processResult(result, outFile,totalStats,classStats);
			if ( isSymmetric)
				withSymm++;

		}
		return withSymm;
	}

	public static CensusResults getExistingResults(String filePath) throws IOException {

		File f = new File(filePath); 

		if (f.exists()) {
			String xml = FileUtils.readFileAsString(filePath);
			CensusResults data = CensusResults.fromXML(xml);

			System.out.println("read " + data.getData().size() + " results from disk...");
			return data;
		}
		return null;
	}


	private boolean processResult(CensusResult result , SynchronizedOutFile outFile, Map<Character, Integer> totalStats,Map<Character,Integer> classStats) throws IOException{

		boolean isSymmetric = false;

		if ( result.getIsSignificant()){
			isSymmetric = true;

		}

		Character scopClass = result.getScopClass();

		String html = ResultConverter.toHTML(result);
		//System.out.println(str.toString());
		outFile.write(html + newline);
		trackStats(totalStats,scopClass,1);
		if ( isSymmetric) {
			trackStats(classStats,scopClass,1);
		}
		return isSymmetric;

	}

	private void trackStats(Map<Character, Integer> totalStats,
			Character scopClass, int i) {

		Integer number = totalStats.get(scopClass);
		if ( number == null) {
			number = 0;

		}

		number += i;
		totalStats.put(scopClass, number);
	}
	private void writeResults(CensusResults census, String filePath){

		String xml = census.toXML();
		//System.out.println(xml);
		FileUtils.writeStringToFile(xml, filePath);
		System.out.println("wrote " + census.getData().size() + " results to disk...");

	}
	/** return the names of the domains that we already analyzed
	 * 
	 * @param census
	 * @return
	 */
	private List<String> getKnownResults(CensusResults census) {

		List<String> known = new ArrayList<String>();
		List<CensusResult> data = census.getData();

		for (CensusResult d : data){
			known.add(d.getName());
		}

		return known;



	}

}
