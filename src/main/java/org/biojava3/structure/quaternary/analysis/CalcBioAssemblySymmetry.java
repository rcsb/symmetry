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
package org.biojava3.structure.quaternary.analysis;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.ResourceManager;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.QuatSymmetryResults;
import org.biojava3.structure.quaternary.core.QuatSymmetryDetector;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.core.Subunits;
import org.biojava3.structure.quaternary.jmolScript.JmolSymmetryScriptGenerator;

public class CalcBioAssemblySymmetry {
	private Structure bioAssembly;
	private QuatSymmetryParameters parameters;
	private QuatSymmetryResults globalSymmetry;
	private AxisTransformation axisTransformation;

	private JmolSymmetryScriptGenerator scriptGenerator;
	
	static public String version;
	static public String build; 
	static {
		try {
			ResourceManager about = ResourceManager.getResourceManager("aboutplayground");

			version = about.getString("project_version");
			build   = about.getString("build");

		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public CalcBioAssemblySymmetry(Structure bioAssembly, QuatSymmetryParameters parameters){
		this.bioAssembly = bioAssembly;
		this.parameters = parameters;
	}
	
	public boolean orient(){
		QuatSymmetryDetector detector = new QuatSymmetryDetector(bioAssembly, parameters);	
		boolean hasProtein = detector.hasProteinSubunits();

		if (hasProtein) {		
			globalSymmetry = detector.getGlobalSymmetry();
			if (parameters.isVerbose()) {
				System.out.println();
				System.out.println("Results for " + Math.round(parameters.getSequenceIdentityThreshold()*100) + "% sequence identity threshold:");
				System.out.println();
				System.out.println("Global symmetry:");
				System.out.println("Stoichiometry         : " + globalSymmetry.getSubunits().getStoichiometry());
				System.out.println("Pseudostoichiometry   : " + globalSymmetry.getSubunits().isPseudoStoichiometric());
				System.out.println("Min sequence identity : " + Math.round(globalSymmetry.getSubunits().getMinSequenceIdentity()*100));
				System.out.println("Max sequence identity : " + Math.round(globalSymmetry.getSubunits().getMaxSequenceIdentity()*100));
				System.out.println("Point group           : " + globalSymmetry.getRotationGroup().getPointGroup());				
				System.out.println("Symmetry RMSD         : " + (float) globalSymmetry.getRotationGroup().getAverageTraceRmsd());
			}

			axisTransformation = new AxisTransformation(globalSymmetry);

			// use factory method to get point group specific instance of script generator
			scriptGenerator = JmolSymmetryScriptGenerator.getInstance(axisTransformation, "g");
			hasProtein = true;
			
		} else {
			System.out.println("No protein chains found for " + bioAssembly.getPDBCode() );
		}
		
		if (parameters.isVerbose() && detector.getLocalSymmetry().size() > 0) {
			System.out.println();
			System.out.println("Local symmetry: ");
			int count = 0;
			for (QuatSymmetryResults localSymmetry: detector.getLocalSymmetry()) {
				AxisTransformation at = new AxisTransformation(localSymmetry);
				System.out.println();
				System.out.println("Results for " + Math.round(parameters.getSequenceIdentityThreshold()*100) + "% sequence identity threshold:");
				System.out.println("Stoichiometry         : " + localSymmetry.getSubunits().getStoichiometry());
				System.out.println("Pseudostoichiometry   : " + localSymmetry.getSubunits().isPseudoStoichiometric());
				System.out.println("Min sequence identity : " + Math.round(localSymmetry.getSubunits().getMinSequenceIdentity()*100));
				System.out.println("Max sequence identity : " + Math.round(localSymmetry.getSubunits().getMaxSequenceIdentity()*100));
				System.out.println("Point group           : " + localSymmetry.getRotationGroup().getPointGroup());				
				System.out.println("Symmetry RMSD         : " + (float) localSymmetry.getRotationGroup().getAverageTraceRmsd());
				System.out.println();
				JmolSymmetryScriptGenerator gen = JmolSymmetryScriptGenerator.getInstance(at, "l"+count);
				if (count == 0) {
					System.out.println(gen.getDefaultOrientation());
				}
				System.out.println(gen.drawPolyhedron());
				System.out.println(gen.drawAxes());
				System.out.println(gen.colorBySymmetry());

				count++;
			}
			System.out.println("draw poly* on; draw axes* on;");
		}
		return hasProtein;
	}
	

	public RotationGroup getRotationGroup() {
		return globalSymmetry.getRotationGroup();
	}

	public Subunits getSubunits() {
		return globalSymmetry.getSubunits();
	}
	
	public AxisTransformation getAxisTransformation() {
		return axisTransformation;
	}

	public JmolSymmetryScriptGenerator getScriptGenerator() {
		return scriptGenerator;
	}
}
