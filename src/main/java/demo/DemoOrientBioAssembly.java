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
import org.biojava3.structure.quaternary.jmolScript.JmolSymmetryScriptGenerator;

public class DemoOrientBioAssembly {
	
	public static void main(String[] args){

		//String[] pdbIDs = new String[]{"4INU", "4D8s","4EAR","4IYQ","3ZKR",};
		String[] pdbIDs = new String[]{"4F88",};
		
		/*
			 2WPD has 2 local symmetries. I’ve attached my presentation from last week.
			 
			 
			Other examples with a single local symmetry are:
			4F88 – local C8
			1LTI – local C5
			2W6E – local C3
			2LXC – local C2
			3OE7 – local C3
			
		*/
		
		for ( String pdbID : pdbIDs){
			runPDB(pdbID);
		}

	}

	public static void runPDB(String pdbID){

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
			
			s = StructureIO.getBiologicalAssembly(pdbID, biolAssemblyNr);

			// Alternative access to structure:			
			//
			//s = readStructure(pdbID, biolAssemblyNr);
			
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
		QuatSymmetryParameters parameters = new QuatSymmetryParameters();
		parameters.setVerbose(false);
		parameters.setSequenceIdentityThreshold(threshold);
		
		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry(s, parameters);

		QuatSymmetryDetector detector = calc.orient();
		
		boolean hasProtein = detector.hasProteinSubunits();

		if (hasProtein) {
			

			String script = "set defaultStructureDSSP true; set measurementUnits ANGSTROMS;  select all;  spacefill off; wireframe off; " +
					"backbone off; cartoon on; color cartoon structure; color structure;  select ligand;wireframe 0.16;spacefill 0.5; " +
					"color cpk ; select all; model 0;set antialiasDisplay true; autobond=false;save STATE state_1;" ;

			String jmolScript = "";
			if ( detector.getLocalSymmetry().size() > 0){
				int count = 0;
				for (QuatSymmetryResults localSymmetry: detector.getLocalSymmetry()) {
					AxisTransformation at = new AxisTransformation(localSymmetry);
					System.out.println();
					System.out.println("Results for " + Math.round(parameters.getSequenceIdentityThreshold()*100) + "% sequence identity threshold:");
					System.out.println("Stoichiometry       : " + localSymmetry.getSubunits().getStoichiometry());
					System.out.println("Pseudostoichiometry : " + localSymmetry.getSubunits().isPseudoStiochiometric());
					System.out.println("Point group         : " + localSymmetry.getRotationGroup().getPointGroup());				
					System.out.println("Symmetry RMSD       : " + (float) localSymmetry.getRotationGroup().getAverageTraceRmsd());
					System.out.println();
					JmolSymmetryScriptGenerator gen = JmolSymmetryScriptGenerator.getInstance(at, "l"+count);
					if (count == 0) {
						script +=  gen.getDefaultOrientation();
					}
					script +=  gen.drawPolyhedron();
					script += gen.drawAxes();
					script +=  gen.colorBySymmetry();

					count++;
				}
				script += "draw poly* on; draw axes* on;";
			} else {
				jmolScript = calc.getScriptGenerator().playOrientations();
			}
			
			
			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(s);

			String title = "Symmetry results for " + pdbID + " bio assembly: " + biolAssemblyNr + " seq cutoff:" + parameters.getSequenceIdentityThreshold();
			jmol.setTitle(title);
			jmol.evalString(script);
			jmol.evalString(jmolScript);
			
			
			
			
			// items for database:
			/*
			System.out.println("=================");
			System.out.println(title );
			System.out.println("=================");
			System.out.println("Sequence ID        : " + parameters.getSequenceIdentityThreshold() );
			System.out.println("Stoichiometry      : " + calc.getSubunits().getStoichiometry());
			System.out.println("Pseudostoichiometry: " + calc.getSubunits().isPseudoStiochiometric());
			System.out.println("Point Group        : " + calc.getRotationGroup().getPointGroup());

			System.out.println("Symmetry RMSD      : " + String.format("%.2f",calc.getRotationGroup().getAverageTraceRmsd()));

			System.out.println("Transf. matrix     : " + calc.getAxisTransformation().getTransformation());
			System.out.println("Geomet. tansf      : " + calc.getAxisTransformation().getGeometicCenterTransformation());
			System.out.println("Dimension          : " + calc.getAxisTransformation().getDimension());
			System.out.println("Subunit count      : " + calc.getSubunits().getSubunitCount());

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
	*/
			System.out.println("=================");
			return calc.getSubunits().isPseudoStiochiometric();
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
