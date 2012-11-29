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

		String pdbID = "4hhb";
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
		
		String jmolScript = calc.animate();

		String script = "set defaultStructureDSSP true; set measurementUnits ANGSTROMS;  select all;  spacefill off; wireframe off; " +
				"backbone off; cartoon on; color cartoon structure; color structure;  select ligand;wireframe 0.16;spacefill 0.5; " +
				"color cpk ; select all; model 0;set antialiasDisplay true; autobond=false;save STATE state_1;" ;


		StructureAlignmentJmol jmol = new StructureAlignmentJmol();
		jmol.setStructure(s);

		String title = "Symmetry results for " + pdbID + " bio assembly: " + biolAssemblyNr + " seq cutoff:" + calc.getParams().getSequenceIdentityThreshold();
		jmol.setTitle(title);
		
		
		jmol.evalString(script);


		jmol.evalString(jmolScript);

		
		
		System.out.println("=================");
		System.out.println(title );
		System.out.println("=================");
		System.out.println("Sequence ID   : " + calc.getParams().getSequenceIdentityThreshold() );
		System.out.println("Stoichiometry : " + calc.getFinder().getCompositionFormula());
		System.out.println("Point Group   : " + calc.getRotationGroup().getPointGroup()	);
		System.out.println("Symmetry RMSD : " + String.format("%.2f",calc.getRotationGroup().getAverageTraceRmsd()));
	
		System.out.println("transf. matrix: " + calc.getAxisTransformation().getTransformation());
		
		System.out.println("dimension     : " + calc.getAxisTransformation().getDimension());
		
		System.out.println("=================");
		
	}
	
}
