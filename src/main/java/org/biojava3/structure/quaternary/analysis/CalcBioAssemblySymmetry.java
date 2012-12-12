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

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.FindQuarternarySymmetry;
import org.biojava3.structure.quaternary.core.QuatSymmetryParameters;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.jmolScript.JmolSymmetryScriptGenerator;

public class CalcBioAssemblySymmetry {
	Structure bioAssembly;
	QuatSymmetryParameters params = new QuatSymmetryParameters();
	
	FindQuarternarySymmetry finder;
	RotationGroup rotationGroup;
	AxisTransformation axisTransformation;
	private JmolSymmetryScriptGenerator scriptGenerator;
	
	
	public static String version = "0.0.6";
	
	public CalcBioAssemblySymmetry(){
	}
	
	
	public boolean orient(){
		boolean hasProtein = false;
		finder = new FindQuarternarySymmetry(bioAssembly, params);	

		rotationGroup = new RotationGroup();
		if (finder.getChainCount() > 0) {
			System.out.println();
			rotationGroup = finder.getRotationGroup();
			
			if (params.isVerbose()) {
				System.out.println("Results for " + Math.round(params.getSequenceIdentityThreshold()*100) + "% sequence identity threshold:");
				System.out.println("Stoichiometry  : " + finder.getCompositionFormula());
				System.out.println("Point group    : " + rotationGroup.getPointGroup());	
				System.out.println("Pseudosymmetric: " + finder.isPseudoSymmetric());	
				System.out.println("Symmetry RMSD  : " + (float) rotationGroup.getAverageTraceRmsd());
			}

			axisTransformation = new AxisTransformation(finder.getSubunits(), rotationGroup);

			// use factory method to get point group specific instance of script generator
			scriptGenerator = JmolSymmetryScriptGenerator.getInstance(axisTransformation);
			hasProtein = true;
		}
		return hasProtein;
	}
	
	public Structure getBioAssembly() {
		return bioAssembly;
	}


	public void setBioAssembly(Structure bioAssembly) {
		this.bioAssembly = bioAssembly;
	}


	public QuatSymmetryParameters getParams() {
		return params;
	}


	public void setParams(QuatSymmetryParameters params) {
		this.params = params;
	}


	public FindQuarternarySymmetry getFinder() {
		return finder;
	}


	public void setFinder(FindQuarternarySymmetry finder) {
		this.finder = finder;
	}


	public RotationGroup getRotationGroup() {
		return rotationGroup;
	}


	public void setRotationGroup(RotationGroup rotationGroup) {
		this.rotationGroup = rotationGroup;
	}


	public AxisTransformation getAxisTransformation() {
		return axisTransformation;
	}


	public void setAxisTransformation(AxisTransformation axistTransformation) {
		this.axisTransformation = axistTransformation;
	}


	public JmolSymmetryScriptGenerator getScriptGenerator() {
		return scriptGenerator;
	}


	public void setScriptGenerator(JmolSymmetryScriptGenerator scriptGenerator) {
		this.scriptGenerator = scriptGenerator;
	}


	public void destroy() {
		// clean up all references so this can get garbage collected
		bioAssembly = null;
		finder = null;
		rotationGroup = null;
		axisTransformation = null;
		
	}
	
	
	
}
