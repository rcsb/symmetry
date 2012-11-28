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
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.quaternary.analysis.CalcBioAssemblySymmetry;

public class DemoOrientBioAssembly {
	public static void main(String[] args){

		String pdbID = "4hhb";
		int biolAssemblyNr = 1;
		
		Structure s;
		try {
			s = StructureIO.getBiologicalAssembly(pdbID, biolAssemblyNr);

			CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry();

			calc.setBioAssembly(s);

			calc.orient();
	
			String jmolScript = calc.animate();
			
			
			
			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(s);
			jmol.evalString(StructureAlignmentJmol.DEFAULT_SCRIPT);
			jmol.evalString("backbone off; cartoon on; color structure");
			
	
			jmol.evalString(jmolScript);

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 



	}
}
