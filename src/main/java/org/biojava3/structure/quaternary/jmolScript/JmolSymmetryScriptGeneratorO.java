/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.geometry.Octahedron;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorO extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorO(AxisTransformation axisTransformation, RotationGroup rotationGroup) {
		super(axisTransformation, rotationGroup);
		Octahedron o = new Octahedron();
		double radius = Math.max(axisTransformation.getDimension().z, axisTransformation.getXYRadius());
		o.setMidRadius(radius);
		polyhedron = o;
	}
	
	public int getDefaultZoom() {
		return 60;
	}
	
}
