/**
 * 
 */
package org.biojava3.structure.align.symm.jmolScript;

import org.biojava3.structure.align.symm.geometry.RectangularPrism;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;

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
