/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.geometry.RectangularPrism;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorC1 extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorC1(AxisTransformation axisTransformation, String name) {
		super(axisTransformation, name);
		setPolyhedron(new RectangularPrism(axisTransformation.getDimension().z*2, axisTransformation.getDimension().x*2, axisTransformation.getDimension().y*2));
	}
	
	public int getZoom() {
		// find maximum extension of structure
		double maxExtension = getMaxExtension();
		// find maximum extension of polyhedron
		AxisTransformation at = getAxisTransformation();
		double polyhedronExtension = Math.max(at.getDimension().x, at.getDimension().y);
		
		polyhedronExtension = Math.max(at.getDimension().z, polyhedronExtension);
		int zoom = Math.round((float)(maxExtension/polyhedronExtension * 110));
		if (zoom > 100) {
			zoom = 100;
		}
		return zoom;
	}
	
	public int getOrientationCount() {
		// the last two views (top, bottom) are not that interesting.
		return getPolyhedron().getViewCount()-2;
	}

}
