/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.geometry.Prism;


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
