/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.geometry.Icosahedron;


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
