/**
 * 
 */
package org.biojava3.structure.align.symm.quaternary;

import java.util.List;

import javax.vecmath.Point3d;

/**
 * @author Peter
 *
 */
public interface Polyhedron {

	public Point3d[] getVertices();
	public List<int[]> getLineLoops();
	
}
