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

import java.util.List;

import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava3.structure.align.symm.census.CensusResult;
import org.biojava3.structure.align.symm.census.CensusResults;

/** Do Database Searches for all Protodomains
 * 
 * @author Andreas Prlic
 *
 */
public class ScanProtodomainsDB {

	public static void main(String[] args){
		String path =  "/Users/ap3/WORK/PDB/";

		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,path);

		ScanProtodomainsDB me = new ScanProtodomainsDB();

		me.run();

	}


	public void run(){
		try {
			CensusResults census = CensusResults.getExistingResults();
			
			List<CensusResult> data = census.getData();
			
			for ( CensusResult cr : data){

				// we only work with significant results
				if ( ! cr.getIsSignificant()) {
					continue;
				}
				
				
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
