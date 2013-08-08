/**
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on Sep 9, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava3.structure.align.symm.census;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Future;

import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.SynchronizedOutFile;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.core.util.ConcurrencyTools;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.utils.FileUtils;


/**
 * @deprecated See {@link Census} instead
 */
@Deprecated
public class ScanSCOPForSymmetry implements Runnable {

	public static String newline = System.getProperty("line.separator");
	//protected static ExecutorService pool;

	static {
		int maxThreads = Runtime.getRuntime().availableProcessors()-1;
		if ( maxThreads < 1)
			maxThreads = 1;

		System.out.println("using " + maxThreads + " threads for alignment calculation of DB jobs...");
		//pool = Executors.newFixedThreadPool(maxThreads);
		ConcurrencyTools.setThreadPoolSize(maxThreads);
	}


	public static void main(String[] args){
		
		try {
			ScanSCOPForSymmetry me = new ScanSCOPForSymmetry();

			ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75B);
			
			me.run();

		} catch (Exception e){
			e.printStackTrace();
		}

	}





	@Override
	public void run() {
		AtomCache cache = new AtomCache();
		
		ScopDatabase scop = ScopFactory.getSCOP();

		List<ScopDescription>superfamilies = scop.getByCategory(ScopCategory.Superfamily);

		System.out.println("got " + superfamilies.size() + " superfamilies");

		String htmlFilename = cache.getPath() + File.separator + "scopCensus.html";

		String xmlFilename = cache.getPath() + File.separator + "scopCensus.xml";
		try {


			CensusResults census = getExistingResults(xmlFilename) ;
			if ( census == null)
				census = new CensusResults();

			List<String> knownResults = getKnownResults(census);
			int count = 0;

			List<Future<CensusResult>> futureData = new ArrayList<Future<CensusResult>>();

			for (ScopDescription superfamily : superfamilies){

				Character scopClass = superfamily.getClassificationId().charAt(0);

				if ( scopClass > 'f')
					continue;

				count++;
				int sunid = superfamily.getSunID();
				List<ScopDomain> familyMembers = scop.getScopDomainsBySunid(sunid);
				
				for (int i = 0; i < 20; i++) {
				
					if (i >= familyMembers.size()) break;
					
				ScopDomain first = familyMembers.get(i);

				String name = first.getScopId();
				if ( knownResults.contains(name))
					continue;

				if ( name.equals("ds046__"))
					continue;

				ScopSymmetryCalculation calc = new ScopSymmetryCalculation();
				calc.setDomain(first);
				calc.setCache(cache);
				calc.setScopDescription(superfamily);
				calc.setCount(count);


				Future<CensusResult> result = ConcurrencyTools.submit(calc);
				futureData.add(result);
				
				}
			}


			File o = new File(htmlFilename);
			if ( o.exists())
				o.delete();
			SynchronizedOutFile htmlOutFile = new SynchronizedOutFile(o);


			ResultConverter.printHTMLHeader(htmlOutFile);

			//int fragmentLength = 8;

			int withSymm = 0;
			Map<Character,Integer> classStats = new HashMap<Character, Integer>();
			Map<Character,Integer> totalStats = new HashMap<Character, Integer>();

			List<CensusResult> allResults = census.getData();

			withSymm = countNrSymm(htmlOutFile, withSymm, classStats, totalStats, allResults);

			System.out.println("starting to process multi threaded results...");
			for (Future<CensusResult> futureResult : futureData) {

				CensusResult result = futureResult.get();
				System.out.println("got result: " + result);
				if ( result == null)
					continue;

				boolean isSymmetric = processResult(result, htmlOutFile,totalStats,classStats);
				if ( isSymmetric)
					withSymm++;

				allResults.add(result);
				//if ( withSymm > 5)
				//	break;
				if ( allResults.size() % 100 == 0)
					writeResults(census, xmlFilename);
			} 

			htmlOutFile.write("</tbody></table>");

			census.setData(allResults);

			writeResults(census, xmlFilename);


			System.out.println("===");
			System.out.println("Overall symmetry: " + String.format("%.2f",(withSymm/(float)count)) + "%" );
			System.out.println("Statistics for SCOP classes:");
			for (Character scopClass: totalStats.keySet()){
				Integer total = totalStats.get(scopClass);
				Integer symm  = classStats.get(scopClass);
				if ( total == null)
					total = 0;
				if ( symm == null)
					symm = 0;
				System.out.println(scopClass);
				System.out.println("Class: " + scopClass + " " + String.format("%.2f",(symm/(float)total))  + "%");
			}

			htmlOutFile.close();
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



	private void writeResults(CensusResults census, String filePath){

		String xml = census.toXML();
		//System.out.println(xml);
		FileUtils.writeStringToFile(xml, filePath);
		System.out.println("wrote " + census.getData().size() + " results to disk...");

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

	public static StructureAlignment getCeSymm(){



		CeSymm ceSymm = new CeSymm();
		StructureAlignmentFactory.addAlgorithm(ceSymm);
		CeParameters params = (CeParameters) ceSymm.getParameters();
		if ( params == null) {
			params = new CeParameters();
			ceSymm.setParameters(params);
		}



		// here how to change the aa subst matrix, SDM is the default matrix
		String matrixName = "PRLA000101";
		//SubstitutionMatrix<AminoAcidCompound> sdm = SubstitutionMatrixHelper.getMatrixFromAAINDEX(matrixName);			
		//params.setSubstitutionMatrix(sdm);
		//SubstitutionMatrix<AminoAcidCompound> max = SubstitutionMatrixHelper.getBlosum85();
		//params.setSubstitutionMatrix(max);		

		// we over-weight sequence
		//params.setSeqWeight(2.0);

		//ceSymm.setParameters(params);

		return ceSymm;
	}
}

