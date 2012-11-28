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

public class CalcBioAssemblySymmetry {
	Structure bioAssembly;
	QuatSymmetryParameters params = OrientBiologicalAssembly.getDefaultParameters();
	
	FindQuarternarySymmetry finder ;
	RotationGroup rotationGroup;
	AxisTransformation axistTransformation;
	
	public CalcBioAssemblySymmetry(){
	}
	
	
	public  void orient(){
		
		 finder = new FindQuarternarySymmetry(bioAssembly, params);	

		 rotationGroup = new RotationGroup();
		if (finder.getChainCount() > 0) {
			System.out.println();
			rotationGroup = finder.getRotationGroup();
			System.out.println("Results for " + Math.round(OrientBiologicalAssembly.SEQUENCE_IDENTITY_THRESHOLD*100) + "% sequence identity threshold:");
			System.out.println("Stoichiometry: " + finder.getCompositionFormula());
			System.out.println("Point group  : " + rotationGroup.getPointGroup());		
			System.out.println("Symmetry RMSD: " + (float) rotationGroup.getAverageTraceRmsd());
		} 
					
		axistTransformation = new AxisTransformation(finder.getSubunits(), rotationGroup);
		
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


	public AxisTransformation getAxistTransformation() {
		return axistTransformation;
	}


	public void setAxistTransformation(AxisTransformation axistTransformation) {
		this.axistTransformation = axistTransformation;
	}


	public void destroy() {
		// clean up all references so this can get garbage collected
		bioAssembly = null;
		finder = null;
		rotationGroup = null;
		axistTransformation = null;
		
	}
	
	
	
}
