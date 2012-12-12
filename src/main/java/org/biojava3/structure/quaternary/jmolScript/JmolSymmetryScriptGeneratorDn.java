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
public class JmolSymmetryScriptGeneratorDn extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorDn(AxisTransformation axisTransformation) {
		super(axisTransformation);
		if (axisTransformation.getRotationGroup().getPointGroup().equals("D2")) {
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
		double polyhedronExtension = Math.max(getPolyhedron().getCirumscribedRadius(), getAxisTransformation().getDimension().z);

		int zoom = Math.round((float)(maxExtension/polyhedronExtension * 110));
		if (zoom > 100) {
			zoom = 100;
		}
		return zoom;
	}
	
}
