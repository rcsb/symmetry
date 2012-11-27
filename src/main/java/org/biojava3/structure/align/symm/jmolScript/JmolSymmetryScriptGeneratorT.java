/**
 * 
 */
package org.biojava3.structure.align.symm.jmolScript;

import org.biojava3.structure.align.symm.geometry.Tetrahedron;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;

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
