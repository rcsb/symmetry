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

public class DemoOrientBioAssembly {
	
	public static void main(String[] args){

		//String[] pdbIDs = new String[]{"4INU", "4D8s","4EAR","4IYQ","3ZKR",};
		String[] pdbIDs = new String[]{"4INU",};
		for ( String pdbID : pdbIDs){
			runPDB(pdbID);
		}

	}

	private static void runPDB(String pdbID){

		pdbID = pdbID.toLowerCase();
		
		int  biolAssemblyNr =1;

		Structure s;
		try {

			//			
			AtomCache cache = new AtomCache();
			FileParsingParameters params = cache.getFileParsingParams();
			params.setAlignSeqRes(true);
			params.setParseCAOnly(false);

			StructureIO.setAtomCache(cache);
			
			//s = StructureIO.getBiologicalAssembly(pdbID, biolAssemblyNr);

			// Alternative access to structure:			
			//
			s = readStructure(pdbID, biolAssemblyNr);
			
			System.out.println("MODELS:" + s.nrModels());
			
			boolean pseudosymmetric = analyzeSymmetry(s,pdbID, biolAssemblyNr, 0.30);

			if (pseudosymmetric) {
				analyzeSymmetry(s,pdbID, biolAssemblyNr, 0.95);
			}
			
		
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}
	
	private static boolean analyzeSymmetry(Structure s,String pdbID, int biolAssemblyNr, double threshold) {

		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry();
		calc.getParams().setVerbose(true);
		calc.setBioAssembly(s);

		calc.getParams().setSequenceIdentityThreshold(threshold);

		boolean hasProtein = calc.orient();

		if (hasProtein) {
			String jmolScript = calc.getScriptGenerator().playOrientations();

			String script = "set defaultStructureDSSP true; set measurementUnits ANGSTROMS;  select all;  spacefill off; wireframe off; " +
					"backbone off; cartoon on; color cartoon structure; color structure;  select ligand;wireframe 0.16;spacefill 0.5; " +
					"color cpk ; select all; model 0;set antialiasDisplay true; autobond=false;save STATE state_1;" ;


			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(s);

			String title = "Symmetry results for " + pdbID + " bio assembly: " + biolAssemblyNr + " seq cutoff:" + calc.getParams().getSequenceIdentityThreshold();
			jmol.setTitle(title);


			jmol.evalString(script);


			jmol.evalString(jmolScript);


			// items for database:

			System.out.println("=================");
			System.out.println(title );
			System.out.println("=================");
			System.out.println("Sequence ID   : " + calc.getParams().getSequenceIdentityThreshold() );
			System.out.println("Stoichiometry : " + calc.getFinder().getCompositionFormula());
			System.out.println("Point Group   : " + calc.getRotationGroup().getPointGroup()	);
			System.out.println("Pseudosymmetry: " + calc.getFinder().isPseudoSymmetric());
			System.out.println("Symmetry RMSD : " + String.format("%.2f",calc.getRotationGroup().getAverageTraceRmsd()));

			System.out.println("Transf. matrix: " + calc.getAxisTransformation().getTransformation());

			System.out.println("Dimension                : " + calc.getAxisTransformation().getDimension());
			System.out.println("Subunit count            : " + calc.getFinder().getChainCount());

			System.out.println("Color by subunit         : " + calc.getScriptGenerator().colorBySubunit());
			System.out.println("Color by subunit length  : " + calc.getScriptGenerator().colorBySubunit().length());
			System.out.println("Color by sequence cluster: " + calc.getScriptGenerator().colorBySequenceCluster());
			System.out.println("Color by seq. clst. len  : " + calc.getScriptGenerator().colorBySequenceCluster().length());
			System.out.println("Color by symmetry        : " + calc.getScriptGenerator().colorBySymmetry());
			System.out.println("Color by symmetry length : " + calc.getScriptGenerator().colorBySymmetry().length());

			System.out.println("Draw axes                : " + calc.getScriptGenerator().drawAxes());
			System.out.println("Draw axes length         : " + calc.getScriptGenerator().drawAxes().length());
			System.out.println("Draw polyhedron          : " + calc.getScriptGenerator().drawPolyhedron());
			System.out.println("Draw polyhedron length   : " + calc.getScriptGenerator().drawPolyhedron().length());

			System.out.println("Zoom                     : " + calc.getScriptGenerator().getZoom());
			System.out.println("Default orientation      : " + calc.getScriptGenerator().getDefaultOrientation());
			System.out.println("Orientation count        : " + calc.getScriptGenerator().getOrientationCount());
			for (int i = 0; i <  calc.getScriptGenerator().getOrientationCount(); i++) {
				System.out.println("Orientation name " + i + "       : " + calc.getScriptGenerator().getOrientationName(i));
				System.out.println("Orientation " + i + "            : " + calc.getScriptGenerator().getOrientation(i));
			}

			System.out.println("=================");
			return calc.getFinder().isPseudoSymmetric();
		} else {
			System.out.println("No protein chains found");
			return false;
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
