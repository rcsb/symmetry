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
 * Created on Nov 27, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package demo;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.quaternary.analysis.CalcBioAssemblySymmetry;
import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.QuatSymmetryDetector;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.QuatSymmetryResults;

import org.biojava3.structure.quaternary.core.SymmetryType;
import org.biojava3.structure.quaternary.jmolScript.JmolSymmetryScriptGenerator;

public class DemoOrientBioAssembly {

	public static void main(String[] args){


		String[] pdbIDs = new String[]{"4HHB","4AQ5","1LTI","1STP","4F88","2W6E","2LXC","3OE7","4INU","4D8s","4EAR","4IYQ","3ZKR"};

		/*		  
		    Local symmetry
		    
			2WPD has 2 local symmetries. 

			Other examples with a single local symmetry are:
			4F88 – local C8
			1LTI – local C5
			2W6E – local C3
			2LXC – local C2
			3OE7 – local C3
			
			Local Pseudosymmetry, structure only
			
			3ZDY

		 */

		for ( String pdbID : pdbIDs)
		{
			runPDB(pdbID);
		}

	}

	public static void runPDB(String pdbID){}

	private static CalcBioAssemblySymmetry  analyzeSymmetry(
			Structure s,
			String pdbID,
			int biolAssemblyNr, 
			double threshold,
			boolean structureOnly) 
	{
		
		QuatSymmetryParameters parameters = new QuatSymmetryParameters();
		
		parameters.setVerbose(false);
		
		
		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry(s, parameters);

		return calc;
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



