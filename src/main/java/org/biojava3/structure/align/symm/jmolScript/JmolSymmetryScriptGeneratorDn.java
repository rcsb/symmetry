/**
 * 
 */
package org.biojava3.structure.align.symm.jmolScript;

import org.biojava3.structure.align.symm.geometry.Prism;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;

/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorDn extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorDn(AxisTransformation axisTransformation, RotationGroup rotationGroup) {
		super(axisTransformation, rotationGroup);
		Prism p = new Prism(rotationGroup.getRotation(0).getFold());
		p.setHeight(axisTransformation.getDimension().z*2);
		p.setInscribedRadius(axisTransformation.getXYRadius());
		polyhedron = p;
	}
	
	public int getDefaultZoom() {
		return 80;
	}
	
}
