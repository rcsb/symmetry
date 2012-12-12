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
import java.io.StringWriter;
import java.util.List;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileConvert;
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
	
			assertEquals(s1.isNmr(),s2.isNmr());
			
			assertEquals("The create structure do not have the same nr. of NMR models", 
						s1.nrModels(), s2.nrModels());

		
			assertEquals(s1.isBiologicalAssembly(), s2.isBiologicalAssembly());
			
			for ( int nrModel = 0 ; nrModel < s1.nrModels(); nrModel ++){
				
				List<Chain> chains1 = s1.getModel(nrModel);
				List<Chain> chains2 = s2.getModel(nrModel);

				assertEquals("the nr of chains in model " + nrModel + " does not match!" , chains1.size(),chains2.size());

				for ( int c = 0 ; c < chains1.size() ; c++){
					
					Chain c1 = s1.getModel(nrModel).get(c);

					Chain c2 = s2.getModel(nrModel).get(c);

					assertEquals("Amino acid toPDB is not equal!" , toPDB(c1), toPDB(c2));

					List<Group> seqres1 = c1.getSeqResGroups();

					List<Group> seqres2 = c2.getSeqResGroups();

					assertEquals("Nr of seqres groups does not match!", seqres1.size(),seqres2.size());

					List<Group> atom1 = c1.getAtomGroups();

					List<Group> atom2 = c2.getAtomGroups();

					assertEquals("Nr of atoms does not match", atom1.size(), atom2.size());

					for ( int g =0 ; g < atom1.size(); g++){
						Group g1 = atom1.get(g);
						Group g2 = atom2.get(g);
						
						assertEquals( g1.getType(),g2.getType());
						
						if ( g1.getType().equals(GroupType.AMINOACID))
							assertEquals(FileConvert.toPDB(g1), FileConvert.toPDB(g2));
						// the numbering can be off for Hetatom groups!!


					}
				}
			}

			
			String stoich1 = getStoichiometry(s1, pdbId, bioAssemblyNr, 0.3);

			
			String stoich2 = getStoichiometry(s2, pdbId, bioAssemblyNr, 0.3);

			assertEquals("the assigned stoichiometries don't match!" ,stoich2, stoich1);

		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}

	}

	private String toPDB(Chain chain){

		StringWriter w = new StringWriter();
		//String chainID = chain.getChainID();
		//if ( chainID.equals(DEFAULTCHAIN) ) chainID = " ";
		// do for all groups
		int nrGroups = chain.getAtomLength();
		for ( int h=0; h<nrGroups;h++){

			Group g= chain.getAtomGroup(h);

			if (g.getType().equals(GroupType.AMINOACID) ){

				String pdb = FileConvert.toPDB(g);
				w.append(pdb);
			}
		}

		return w.toString();
	}



	private static String getStoichiometry(Structure s,String pdbID, int biolAssemblyNr, double threshold) {

		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry();
		calc.getParams().setVerbose(false);
		calc.setBioAssembly(s);

		calc.getParams().setSequenceIdentityThreshold(threshold);

		boolean hasProtein = calc.orient();

		return calc.getFinder().getCompositionFormula();
	}

	private Structure getStructureFromIO(String pdbID, int biolAssemblyNr) throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		params.setAlignSeqRes(true);
		params.setParseCAOnly(false);
		params.setStoreEmptySeqRes(true);
		params.setLoadChemCompInfo(true);
		params.setAtomCaThreshold(Integer.MAX_VALUE);
		cache.setAutoFetch(true);
		params.setParseBioAssembly(true);


		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getBiologicalAssembly(pdbID, biolAssemblyNr);
		return s;
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
