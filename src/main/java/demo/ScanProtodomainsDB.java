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
 * Created on Oct 5, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package demo;

import java.io.File;
import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.gui.AlignmentCalcDB;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava3.core.util.ConcurrencyTools;
import org.biojava3.structure.align.symm.census.CensusResult;
import org.biojava3.structure.align.symm.census.CensusResults;
import org.biojava3.structure.utils.FileUtils;
import org.jmol.adapter.smarter.Atom;
import org.rcsb.fatcat.server.util.ResourceManager;

/** Do Database Searches for all Protodomains
 * 
 * @author Andreas Prlic
 *
 */
public class ScanProtodomainsDB {

	public static void main(String[] args){
		String path =  ResourceManager.getString("pdbFilePath");
		System.out.println("PDB file path:" +  path);

		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,path);

		ScanProtodomainsDB me = new ScanProtodomainsDB();

		me.run();

	}

	
	public static final String RESULTS_DIR="PROTODBSEARCH";
			

	public void run(){
		try {
			CensusResults census = CensusResults.getExistingResults();
			
			List<CensusResult> data = census.getData();
			
			UserConfiguration config = new UserConfiguration();
			AtomCache cache = new AtomCache(config);
			
			String resultsDirS =  config.getPdbFilePath() + File.separator + RESULTS_DIR ;
			
			File resultsDir = new File(resultsDirS);
			if (! resultsDir.exists()){
				
				System.out.println("creating dir " + resultsDir.getAbsolutePath());
				resultsDir.mkdir();
			}
			
			
					
			for ( CensusResult cr : data){

				// we only work with significant results
				if ( ! cr.getIsSignificant()) {
					continue;
				}
				
				String protoDomain = cr.getProtoDomain();
				
				String resultsFileFinalS = resultsDirS + File.separator + protoDomain + File.separator ;
				String fileName =  "results_" + protoDomain + ".out";
			
				File finalFile = new File(resultsFileFinalS + File.separator + fileName);
				
				if ( finalFile.exists() ) {
					System.out.println("not re-running alignments for " + protoDomain);
					continue;
				}
				
				Structure s1 = cache.getStructure(protoDomain);
				
				// create a tmp dir
				
				String outFile = resultsDirS + File.separator  + "tmp";
				File tmpDir = new File(outFile);
				if ( ! tmpDir.exists()){
					tmpDir.mkdir();
				}
				//
				AlignmentCalcDB dbScan = new AlignmentCalcDB(null, s1, protoDomain, config, outFile, config.isSplit());
				
				StructureAlignment algo = getAlgorithm();
				
				dbScan.setAlgorithm(algo);
				int maxSize = ConcurrencyTools.getThreadPool().getMaximumPoolSize();
				dbScan.setNrCPUs(maxSize);
				
				// only one Thread here..
				// dbScan is running in multi threaded mode anyways.
				
				Thread t = new Thread(dbScan);
				t.start();
				
				try {
					Thread.sleep(30000);
					
					
					int queueSize = ConcurrencyTools.getThreadPool().getQueue().size(); 
					// wait until all jobs have finished..
					while ( queueSize > 0  ) {
						
						queueSize = ConcurrencyTools.getThreadPool().getQueue().size(); 
						System.out.println("queue size : " + queueSize + " active: " + ConcurrencyTools.getThreadPool().getActiveCount() + " pool size: " + ConcurrencyTools.getThreadPool().getPoolSize());
						Thread.sleep(1000);

					}
					
				}
				catch (Exception e){
					e.printStackTrace();
					ConcurrencyTools.shutdown();
				}

				//copy results file to Permanent location..
				
				File finalDir = new File(resultsDirS + File.separator + protoDomain);
				if ( ! finalDir.exists()){
					System.out.println("mkdir " + finalDir.getAbsolutePath());
					finalDir.mkdir();
				}
				
				File tmpResult = new File(outFile+File.separator + fileName);
				System.out.println("cp " + tmpResult.getAbsolutePath() + " " + finalFile.getAbsolutePath());
				FileUtils.copy(tmpResult,finalFile);
				
				
				
				
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}


	// using Ce with Seq scoring... 
	private StructureAlignment getAlgorithm() {
		StructureAlignment algo = new CeMain();
		
		CeParameters parameters = new CeParameters();
		parameters.setScoringStrategy(CeParameters.SEQUENCE_CONSERVATION);
		
		return algo;
	}
}
