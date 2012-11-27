/**
 * 
 */
package org.biojava3.structure.align.symm.jmolScript;

import org.biojava3.structure.align.symm.geometry.Icosahedron;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;

/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorI extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorI(AxisTransformation axisTransformation, RotationGroup rotationGroup) {
		super(axisTransformation, rotationGroup);
		double radius = Math.max(axisTransformation.getDimension().z, axisTransformation.getXYRadius());
		Icosahedron i = new Icosahedron();
		i.setMidRadius(radius);
		polyhedron = i;
	}
	
	public int getDefaultZoom() {
		return 80;
	}
	
}
