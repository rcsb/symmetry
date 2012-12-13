/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import org.biojava3.structure.quaternary.core.AxisTransformation;
import org.biojava3.structure.quaternary.geometry.Prism;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorDn extends JmolSymmetryScriptGenerator {

	public JmolSymmetryScriptGeneratorDn(AxisTransformation axisTransformation) {
		super(axisTransformation);
		int fold = axisTransformation.getRotationGroup().getRotation(0).getFold();
		
		// special case for D2. Since there is no 2-fold prism, draw a 4-fold
		// prism that encases the D2 structure
		if (axisTransformation.getRotationGroup().getPointGroup().equals("D2")) {
			fold = 4;
		}
		
		Prism p = new Prism(fold);
		p.setHeight(axisTransformation.getDimension().z*2);
		p.setInscribedRadius(axisTransformation.getXYRadius());
		setPolyhedron(p);
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
