/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.geometry.Prism;
import org.biojava3.structure.quaternary.geometry.RectangularPrism;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorCn extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorCn(AxisTransformation axisTransformation) {
		super(axisTransformation);
		if (axisTransformation.getRotationGroup().getPointGroup().equals("C2")) {
			setPolyhedron(new RectangularPrism(axisTransformation.getDimension().z*2, axisTransformation.getDimension().x*2, axisTransformation.getDimension().y*2));
		} else {
			Prism p = new Prism(axisTransformation.getRotationGroup().getRotation(0).getFold());
			p.setHeight(axisTransformation.getDimension().z*2);
			p.setInscribedRadius(axisTransformation.getXYRadius());
			setPolyhedron(p);
		}
	}
	
	public int getZoom() {
		// find maximum extension of structure
		double maxExtension = getMaxExtension();
		// find maximum extension of polyhedron
		AxisTransformation at = getAxisTransformation();
		double polyhedronExtension = Math.max(getPolyhedron().getCirumscribedRadius(), at.getDimension().z);
		
		int zoom = Math.round((float)(maxExtension/polyhedronExtension * 110));
		if (zoom > 100) {
			zoom = 100;
		}
		return zoom;
	}
	
	public int getOrientationCount() {
		//  the last two views (top, bottom) are not that interesting.
		return getPolyhedron().getViewCount()-2;
	}
}
