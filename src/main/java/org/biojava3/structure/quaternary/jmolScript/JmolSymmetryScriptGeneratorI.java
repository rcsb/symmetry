/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.RotationAxisAligner;
import org.biojava3.structure.quaternary.geometry.Icosahedron;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorI extends JmolSymmetryScriptGeneratorPointGroup {

	public JmolSymmetryScriptGeneratorI(RotationAxisAligner axisTransformation, String name) {
		super(axisTransformation, name);
		double radius = Math.max(axisTransformation.getDimension().z, axisTransformation.getXYRadius());
		Icosahedron i = new Icosahedron();
		i.setMidRadius(radius);
		setPolyhedron(i);
	}
	
	public int getZoom() {
		// find maximum extension of structure
		double maxExtension = getMaxExtension();
		// find maximum extension of polyhedron
		double polyhedronExtension = getPolyhedron().getCirumscribedRadius();
		
		int zoom = Math.round((float)(maxExtension/polyhedronExtension * 110));
		if (zoom > 100) {
			zoom = 100;
		}
		return zoom;
	}
	
}
