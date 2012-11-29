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
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.quaternary.analysis.CalcBioAssemblySymmetry;

public class DemoOrientBioAssembly {
	public static void main(String[] args){

		String pdbID = "1bhc";
		int  biolAssemblyNr = 1;

		Structure s;
		try {

			//			
			AtomCache cache = new AtomCache();
			FileParsingParameters params = cache.getFileParsingParams();
			params.setAlignSeqRes(false);
			params.setParseCAOnly(false);

			StructureIO.setAtomCache(cache);

			s = StructureIO.getBiologicalAssembly(pdbID, biolAssemblyNr);

			analyzeSymmetry(s,pdbID, biolAssemblyNr, 0.30);

			analyzeSymmetry(s,pdbID, biolAssemblyNr, 0.95);
			
		
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 



	}

	private static void analyzeSymmetry(Structure s,String pdbID, int biolAssemblyNr, double threshold) {

		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry();

		calc.setBioAssembly(s);

		calc.getParams().setSequenceIdentityThreshold(threshold);
		
		calc.orient();
		
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
		
		System.out.println("Color by subunit         : " + calc.getScriptGenerator().colorBySubunit());
		System.out.println("Color by sequence cluster: " + calc.getScriptGenerator().colorBySequenceCluster());
		System.out.println("Color by symmetry        : " + calc.getScriptGenerator().colorBySymmetry());
		
		System.out.println("Draw axes                : " + calc.getScriptGenerator().drawAxes());
		System.out.println("Draw polyhedron          : " + calc.getScriptGenerator().drawPolyhedron());
		
		System.out.println("Default orientation  : " + calc.getScriptGenerator().setDefaultOrientation());
		System.out.println("Orientation count        : " + calc.getScriptGenerator().getOrientationCount());
		for (int i = 0; i <  calc.getScriptGenerator().getOrientationCount(); i++) {
			System.out.println("Orientation name " + i + "       : " + calc.getScriptGenerator().getOrientationName(i));
			System.out.println("Orientation " + i + "            : " + calc.getScriptGenerator().setOrientation(i));
		}
		
		System.out.println("=================");
		
	}
	
}
