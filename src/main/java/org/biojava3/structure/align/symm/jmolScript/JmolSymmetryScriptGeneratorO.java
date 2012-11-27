/**
 * 
 */
package org.biojava3.structure.align.symm.jmolScript;

import org.biojava3.structure.align.symm.geometry.Octahedron;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;

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
