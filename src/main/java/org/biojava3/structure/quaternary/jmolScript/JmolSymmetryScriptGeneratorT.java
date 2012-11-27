/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.core.RotationGroup;
import org.biojava3.structure.quaternary.geometry.Tetrahedron;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorT extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorT(AxisTransformation axisTransformation, RotationGroup rotationGroup) {
		super(axisTransformation, rotationGroup);
		double radius = Math.max(axisTransformation.getDimension().z, axisTransformation.getXYRadius());
		Tetrahedron t = new Tetrahedron();
		t.setMidRadius(radius);
		polyhedron = t;
//		System.out.println("Tetrahedron: mid radius" + radius);
//		System.out.println("Tetrahedron: circ. radius" + t.getCirumscribedRadius());
	}
	
	public int getDefaultZoom() {
		return 60;
	}
	
}
