/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.geometry.RectangularPrism;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorC1 extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorC1(AxisTransformation axisTransformation, RotationGroup rotationGroup) {
		super(axisTransformation, rotationGroup);
		polyhedron = new RectangularPrism(axisTransformation.getDimension().z*2, axisTransformation.getDimension().x*2, axisTransformation.getDimension().y*2);
	}
	
	public int getDefaultZoom() {
		return 80;
	}
	
}
