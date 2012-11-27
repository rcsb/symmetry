/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.geometry.Prism;
import org.biojava3.structure.quaternary.geometry.RectangularPrism;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorCn extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorCn(AxisTransformation axisTransformation, RotationGroup rotationGroup) {
		super(axisTransformation, rotationGroup);
		if (rotationGroup.getPointGroup().equals("C2")) {
			polyhedron = new RectangularPrism(axisTransformation.getDimension().z*2, axisTransformation.getDimension().x*2, axisTransformation.getDimension().y*2);
		} else {
			Prism p = new Prism(rotationGroup.getRotation(0).getFold());
			p.setHeight(axisTransformation.getDimension().z*2);
			p.setInscribedRadius(axisTransformation.getXYRadius());
			polyhedron = p;
		}
	}
	
	public int getDefaultZoom() {
		return 80;
	}
	
}
