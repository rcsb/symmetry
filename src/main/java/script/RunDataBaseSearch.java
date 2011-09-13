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
 * Created on Sep 8, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package script;

import java.io.File;
import java.util.SortedSet;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;

import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.align.CallableStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.client.PdbPair;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.SynchronizedOutFile;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.core.util.ConcurrencyTools;
import org.biojava3.structure.align.symm.CeSymm;
import org.rcsb.fatcat.server.PdbChainKey;
import org.rcsb.fatcat.server.dao.SequenceClusterDAO;

import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

/** a Utility method to run a multi-threaded database search
 * 
 * @author Andreas Prlic
 *
 */
public class RunDataBaseSearch implements Runnable{

	
	UserConfiguration config;
	PdbChainKey query ;
	StructureAlignment algorithm;
	int nrCPUs;
	File outPutDir;

	public static final int DEFAULT_CLUSTER_CUTOFF = 40;


	
	public RunDataBaseSearch(PdbChainKey key, StructureAlignment algorithm) {
		query = key;
		this.algorithm = algorithm;
		config = new UserConfiguration();
		nrCPUs = Runtime.getRuntime().availableProcessors()-1;

		if (nrCPUs < 1)
			nrCPUs = 1;

		outPutDir = new File(config.getPdbFilePath() + "dbsearch_" + query.toName() + "_" + algorithm.getAlgorithmName()) ;
		if ( ! outPutDir.exists());
		outPutDir.mkdir();

	}


	public void run() {

		SequenceClusterDAO dao = new SequenceClusterDAO();

		SortedSet<PdbChainKey> representatives = dao.getClusterEntities(DEFAULT_CLUSTER_CUTOFF);

		AtomCache cache = new AtomCache(config);

		cache.getFileParsingParams().setUpdateRemediatedFiles(true);
		
		String outFile = outPutDir+File.separator+query.toName()+".results.out";
		
		System.out.println("writing results to " + outFile);
		File outFileF = new File(outFile);


		
		
		ConcurrencyTools.setThreadPoolSize(nrCPUs);
		
		ThreadPoolExecutor  pool = ConcurrencyTools.getThreadPool();
		
		

		long startTime = System.currentTimeMillis();
		
		try {
			SynchronizedOutFile outF = new SynchronizedOutFile(outFileF);
			Atom[] ca1 = cache.getAtoms( query.toName() );

			
			for (PdbChainKey representative : representatives){

				CallableStructureAlignment ali = new CallableStructureAlignment();

				
				PdbPair pair = new PdbPair(query.toName(), representative.toName());
				try {
					ali.setCa1(ca1);
				} catch (Exception e){
					e.printStackTrace();
					ConcurrencyTools.shutdown();
					return;
				}
				
				ali.setCache(cache);
				ali.setAlgorithmName(algorithm.getAlgorithmName());
				ali.setParameters(algorithm.getParameters());
				ali.setPair(pair);
				ali.setOutFile(outF);			
				ali.setOutputDir(outPutDir);
				ali.setCache(cache);

				Future<AFPChain> afpChain = ConcurrencyTools.submit(ali);

			}
			
			try {
				while ( pool.getCompletedTaskCount() < representatives.size()-1  ) {
					//long now = System.currentTimeMillis();
					//System.out.println( pool.getCompletedTaskCount() + " " + (now-startTime)/1000 + " " + pool.getPoolSize() + " " + pool.getActiveCount()  + " " + pool.getTaskCount()  );
//					if ((now-startTime)/1000 > 60) {
//						
//						interrupt();
//						System.out.println("completed: " + pool.getCompletedTaskCount());
//					}

					
					Thread.sleep(1000);
					long timeN = System.currentTimeMillis();
					long total = timeN -startTime;
					double avTime = 0;
					if ( pool.getCompletedTaskCount() > 0 )
						avTime = (total / 1000f / pool.getCompletedTaskCount());
					
					System.out.println("ran " + pool.getCompletedTaskCount() + "/" + pool.getTaskCount() + " alignments. ("+String.format("%.2f", avTime)+" / per alig.) Using " +pool.getPoolSize() + "CPUs." );
				}
				outF.close();
			}
			catch (Exception e){
				e.printStackTrace();
				
			}
			
		} catch (Exception e){
			e.printStackTrace();
		}
		long now = System.currentTimeMillis();
		System.out.println("Calculation took : " + (now-startTime)/1000 + " sec.");

	}





	public UserConfiguration getUserConfiguration() {
		return config;
	}


	public void setUserConfiguration(UserConfiguration config) {
		this.config = config;
	}


	public int getNrCPUs() {
		return nrCPUs;
	}


	public void setNrCPUs(int nrCPUs) {
		this.nrCPUs = nrCPUs;
	}


	

	
	public static void main(String[] args){

		String name1 = "1hiv.A";

		PdbChainKey key = PdbChainKey.fromName(name1);


		try {
			StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			//StructureAlignment algorithm = getCeSymm();

			RunDataBaseSearch me = new RunDataBaseSearch(key, algorithm);
			
			UserConfiguration config = new UserConfiguration();			
			me.setUserConfiguration(config);
			
			me.run();

		} catch (Exception e){
			e.printStackTrace();
		}

	}
	


}
