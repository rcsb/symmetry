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


import java.io.IOException;
import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
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


		//String[] pdbIDs = new String[]{"4HHB","4AQ5","1LTI","1STP","4F88","2W6E","2LXC","3OE7","4INU","4D8s","4EAR","4IYQ","3ZKR"};
		
		String[] pdbIDs = new String[]{"1LTI"};

		int bioAssemblyNr = 1;
		
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
			try {
			
				runPDB(pdbID,bioAssemblyNr);
			
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
		}

	}

	public static void runPDB(String pdbID, int bioAssemblyNr) throws IOException, StructureException{
		
		
		
		pdbID = pdbID.toLowerCase();
		
		
		//Structure s = StructureIO.getBiologicalAssembly(pdbID, bioAssemblyNr);
		Structure s = readStructure(pdbID, bioAssemblyNr);
		
		QuatSymmetryParameters parameters = new QuatSymmetryParameters();

		parameters.setVerbose(true);


		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry(s, parameters);
		
		QuatSymmetryDetector detector = calc.orient();
		
		List<QuatSymmetryResults> globalResults = detector.getGlobalSymmetry();
		
		System.out.println("# of global results: " + globalResults.size());
		
		List<List<QuatSymmetryResults>> localResults = detector.getLocalSymmetries();
		
		
		
		showResults(s, pdbID + "[" + bioAssemblyNr + "] Global", globalResults);
		
		
		for (int counter = 0;counter < localResults.size() ; counter++){
			List<QuatSymmetryResults> localResultsL = localResults.get(counter);
			
			showResults(s,pdbID + "[" + bioAssemblyNr + "] Local #" + (counter+1) , localResultsL);
		}
		
		
		
	}

	private static void showResults(Structure s, String title,
			List<QuatSymmetryResults> results) {
		
		
		int count = 0 ;
		for (QuatSymmetryResults result: results) {
			
			String script = "set defaultStructureDSSP true; set measurementUnits ANGSTROMS;  select all;  spacefill off; wireframe off; " +
					"backbone off; cartoon on; color cartoon structure; color structure;  select ligand;wireframe 0.16;spacefill 0.5; " +
					"color cpk ; select all; model 0;set antialiasDisplay true; autobond=false;save STATE state_1;" ;
			count++;
			
			if ( result.getSubunits().isPseudoSymmetric()) {
				System.out.println("pseudosymmetric!");
			} else {
				System.out.println(" not pseudosymmetric!");
				
			}
			
			AxisTransformation axisTransformation = new AxisTransformation(result);

			// use factory method to get point group specific instance of script generator
			JmolSymmetryScriptGenerator scriptGenerator = JmolSymmetryScriptGenerator.getInstance(axisTransformation, "g");

			script += scriptGenerator.getOrientationWithZoom(0);
			script += scriptGenerator.drawPolyhedron();
			script += scriptGenerator.drawAxes();
			script += scriptGenerator.colorBySymmetry();
			
			String longTitle = title + " count:"+ count + " [" + result.getSubunits().getStoichiometry() +"]";
			
			script += "draw axes* on; draw poly* on;";
			
			
			// show in Jmol...

			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(s);
			
			jmol.setTitle(longTitle);
			jmol.evalString(script);
		}
		
		
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



