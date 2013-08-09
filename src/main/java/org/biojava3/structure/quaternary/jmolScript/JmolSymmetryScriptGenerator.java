package org.biojava3.structure.quaternary.jmolScript;

import javax.vecmath.Matrix4d;

import org.biojava3.structure.quaternary.core.AxisAligner;
import org.biojava3.structure.quaternary.core.HelixAxisAligner;
import org.biojava3.structure.quaternary.core.RotationAxisAligner;

public abstract class JmolSymmetryScriptGenerator {

	/**
	 * Returns an instance of a JmolSymmetryScriptGenerator, based on the symmetry of a structure (factory method)
	 * @param axisAligner
	 * @param rotationGroup
	 * @return instance of JmolSymmetryScriptGenerator
	 */
	public static JmolSymmetryScriptGenerator getInstance(AxisAligner axisAligner, String name) {
		String symmetry = axisAligner.getSymmetry();
		
		if (symmetry.equals("C1")) {
			return new JmolSymmetryScriptGeneratorC1((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.startsWith("C")) {
			return new JmolSymmetryScriptGeneratorCn((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.startsWith("D")) {
			return new JmolSymmetryScriptGeneratorDn((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.equals("T")) {
			return new JmolSymmetryScriptGeneratorT((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.equals("O")) {
			return new JmolSymmetryScriptGeneratorO((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.equals("I")) {
			return new JmolSymmetryScriptGeneratorI((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.equals("H")) {
			return new JmolSymmetryScriptGeneratorH((HelixAxisAligner)axisAligner, name);
		}
		
		return null;
	}
	/**
	 * Returns the Jmol zoom to fit polyhedron and symmetry axes. This zoom
	 * level should be used so that the polyhedron and symmetry axes are not cutoff.
	 * @return
	 */
	abstract public int getZoom();

	/**
	 * Returns a Jmol script to set the default orientation for a structure
	 * @return Jmol script
	 */
	public abstract String getDefaultOrientation();

	/**
	 * Returns the number of orientations available for this structure
	 * @return number of orientations
	 */
	public abstract int getOrientationCount();

	/**
	 * Returns a Jmol script that sets a specific orientation
	 * @param index orientation index
	 * @return Jmol script
	 */
	public abstract String getOrientation(int index);

	/**
	 * Returns a Jmol script that sets a specific orientation and zoom
	 * to draw either axes or polyhedron
	 * @param index orientation index
	 * @return Jmol script
	 */
	public abstract String getOrientationWithZoom(int index);

	/**
	 * Returns the name of a specific orientation
	 * @param index orientation index
	 * @return name of orientation
	 */
	public abstract String getOrientationName(int index);

	/**
	 * Returns transformation matrix to orient structure
	 * @return transformation matrix
	 */
	public abstract Matrix4d getTransformation();
	
	/** Sets a default Jmol script used for coloring. This method is
	 * used in local symmetry cases to color those subunits that are
	 * not related by symmetry.
	 * @param colorScript
	 */	
	public abstract void setDefaultColoring(String colorScript);
	
	/**
	 * Returns a Jmol script that draws an invisible polyhedron around a structure.
	 * Use showPolyhedron() and hidePolyhedron() to toggle visibility.
	 * @return Jmol script
	 */
	public abstract String drawPolyhedron();

	public abstract String hidePolyhedron();

	public abstract String showPolyhedron();

	/**
	 * Returns a Jmol script that draws symmetry or inertia axes for a structure.
	 * Use showAxes() and hideAxes() to toggle visibility.
	 * @return Jmol script
	 */
	public abstract String drawAxes();

	/**
	 * Returns a Jmol script to hide axes
	 * @return Jmol script
	 */
	public abstract String hideAxes();

	/**
	 * Returns a Jmol script to show axes
	 * @return Jmol script
	 */
	public abstract String showAxes();

	/**
	 * Returns a Jmol script that displays a symmetry polyhedron and symmetry axes
	 * and then loop through different orientations
	 * @return Jmol script
	 */
	public abstract String playOrientations();

	/**
	 * Returns a Jmol script that colors the subunits of a structure by different colors
	 * @return
	 */
	public abstract String colorBySubunit();

	/**
	 * Returns a Jmol script that colors subunits by their sequence cluster ids.
	 * @return Jmol script
	 */
	public abstract String colorBySequenceCluster();

	/**
	 * Returns a Jmol script that colors subunits to highlight the symmetry within a structure
	 * @return Jmol script
	 */
	public abstract String colorBySymmetry();

}