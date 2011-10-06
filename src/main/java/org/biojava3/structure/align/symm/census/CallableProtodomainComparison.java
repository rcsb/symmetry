package org.biojava3.structure.align.symm.census;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;

public class CallableProtodomainComparison implements Callable<ProtoDomainDBSearchResults>{

	AtomCache cache ;

	List<CensusResult> data ;

	String protoDomain;
	

	public static String getResultFilePath(AtomCache cache, String protoDomain){
		String dor = cache.getPath() + File.separator + "protodomain" + File.separator;
		String results = dor + protoDomain + ".dbsearch";

		return results;
	}
	
	@Override
	public ProtoDomainDBSearchResults call() throws Exception {


		
		if ( protoDomain == null || protoDomain.equals("")){
			System.err.println("can't run DB search for emtpy protodomain!");
			return null;
		}

		String dor = cache.getPath() + File.separator + "protodomain" + File.separator;
		File d = new File(dor);
		if ( ! d.exists()) {
			System.out.println("creating dir " + dor);
			d.mkdir();
		}

		// check if result file already exist
		String results = getResultFilePath(cache, protoDomain);
		File resultsFile = new File(results);
		
		if ( resultsFile.exists()) {
			// we already have run this calculation, don't re-run!
			System.out.println("not re-processing " + protoDomain);

			ProtoDomainDBSearchResults dbResults =  ProtoDomainDBSearchResults.fromFile(resultsFile);
			return dbResults;
		}

		// this protodomain has not been scanned yet.
		// run DB scan

		
		Atom[] ca1 = null ;
		try {
			ca1 = cache.getAtoms(protoDomain);
		} catch (Exception e){
			System.out.println("problem loading protoDomain:" + protoDomain);
			return null;
		}

		StructureAlignment algorithm = getAlgorithm();

		List<ProtoDomainDBSearchResult> dbSearchResults = new ArrayList<ProtoDomainDBSearchResult>();

		int counter = 0;
		for ( CensusResult result : data ){
			if ( ! result.getIsSignificant()) {
				continue;
			}
			String protoDomain2 = null;
			try {
				counter ++;
				protoDomain2 = result.getProtoDomain();

				Atom[] ca2 =  null;
				try {
					ca2 = cache.getAtoms(protoDomain2);
				} catch (Exception e){
					System.out.println("problem loading " + protoDomain);
					continue;
				}

				//System.out.print("#" + counter + " " +protoDomain + " vs. " + protoDomain2 + " ");

				AFPChain afpChain= algorithm.align(ca1,ca2);

				if ( afpChain == null)
					continue;

				ProtoDomainDBSearchResult pddsr = convertAFPChain(afpChain,result);

				//System.out.print(pddsr);

				dbSearchResults.add(pddsr);

			} catch (Exception e){
				//e.printStackTrace();
				System.err.println("Could not align " + protoDomain + " " + protoDomain2 + " " + e.getMessage());
			}


		}

		ProtoDomainDBSearchResults finalResults = new ProtoDomainDBSearchResults();
		finalResults.setDbSearchResults(dbSearchResults);
		finalResults.setProtoDomain(protoDomain);
		finalResults.writeToFile(resultsFile);

		return finalResults;

	}



	private ProtoDomainDBSearchResult convertAFPChain(AFPChain afpChain,CensusResult census) {


		ProtoDomainDBSearchResult result = new ProtoDomainDBSearchResult();

		result.setName1(afpChain.getName1());
		result.setName2(afpChain.getName2());
		result.setAlgorithmName(CeCPMain.algorithmName);
		result.setAlignScore(new Double(afpChain.getAlignScore()).floatValue());
		result.setCa1Length(afpChain.getCa1Length());
		result.setCa2Length(afpChain.getCa2Length());
		result.setCoverage1(afpChain.getCoverage1());
		result.setCoverage2(afpChain.getCoverage2());
		result.setEqr(afpChain.getNrEQR());
		result.setIdentity(new Double(afpChain.getIdentity()).floatValue());
		result.setSimilarity(new Double(afpChain.getSimilarity()).floatValue());
		result.setProbability(afpChain.getProbability());
		result.setRmsd(new Double(afpChain.getTotalRmsdOpt()).floatValue());
		result.setOriginalScopID2(census.getName());
		result.setClassificationID2(census.getClassificationId());
		if ( afpChain.getBlockNum() > 1){
			result.setIsCP(true);
		} else {
			result.setIsCP(false);
		}

		return result;





	}



	private StructureAlignment getAlgorithm() {

		CeCPMain algo = new CeCPMain();

		return algo;
	}



	public AtomCache getCache() {
		return cache;
	}



	public void setCache(AtomCache cache) {
		this.cache = cache;
	}



	public List<CensusResult> getData() {
		return data;
	}



	public void setData(List<CensusResult> data) {
		this.data = data;
	}



	public String getProtoDomain() {
		return protoDomain;
	}



	public void setProtoDomain(String protoDomain) {
		this.protoDomain = protoDomain;
	}



}
