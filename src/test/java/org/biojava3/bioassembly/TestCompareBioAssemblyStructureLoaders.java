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
 * Created on Dec 12, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava3.bioassembly;

import java.io.IOException;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.quaternary.analysis.CalcBioAssemblySymmetry;

import junit.framework.TestCase;

public class TestCompareBioAssemblyStructureLoaders 
extends TestCase{

	public void test1STP(){

		String pdbId = "1STP";

		int bioAssemblyNr = 1;

		testLoaders(pdbId,bioAssemblyNr);
	}

	private void testLoaders(String pdbId, int bioAssemblyNr) {

		try {
			Structure s1 = getStructureFromIO(pdbId,bioAssemblyNr);
			Structure s2 = readStructure(pdbId, bioAssemblyNr);
			
			String stoich1 = getStoichiometry(s1, pdbId, bioAssemblyNr, 0.3);
			String stoich2 = getStoichiometry(s2, pdbId, bioAssemblyNr, 0.3);
			
			assertEquals(stoich2, stoich1);
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}

	}

	private Structure getStructureFromIO(String pdbID, int biolAssemblyNr) throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		params.setAlignSeqRes(true);
		params.setParseCAOnly(false);

		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getBiologicalAssembly(pdbID, biolAssemblyNr);
		return s;
	}

	private static String getStoichiometry(Structure s,String pdbID, int biolAssemblyNr, double threshold) {

		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry();
		calc.getParams().setVerbose(true);
		calc.setBioAssembly(s);

		calc.getParams().setSequenceIdentityThreshold(threshold);

		boolean hasProtein = calc.orient();

		return calc.getFinder().getCompositionFormula();
	}
	
	private static Structure  readStructure(String pdbId, int bioAssemblyId) {
		// initialize the PDB_DIR env variable
		AtomCache cache = new AtomCache();

		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		p.setLoadChemCompInfo(true);
		p.setAtomCaThreshold(Integer.MAX_VALUE);
		//p.setAcceptedAtomNames(new String[]{" CA "});
		p.setParseBioAssembly(true);



		PDBFileReader pdbreader = new PDBFileReader();
		pdbreader.setPath(cache.getPath());
		pdbreader.setFileParsingParameters(p);
		pdbreader.setAutoFetch(true);
		pdbreader.setBioAssemblyId(bioAssemblyId);
		pdbreader.setBioAssemblyFallback(false);
		Structure structure = null;
		try { 
			structure = pdbreader.getStructureById(pdbId);
			if ( bioAssemblyId > 0 )
				structure.setBiologicalAssembly(true);
			structure.setPDBCode(pdbId);
		} catch (Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
		return structure;
	}
}
