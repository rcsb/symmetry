package org.biojava.nbio.structure.align.symm.axis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.vecmath.Matrix4d;

/**
 * Data Structure that stores all the symmetry axis that
 * describe the symmetry of a structure. Generalizes to 
 * all types of symmetry, the classic ones (Cn, Dn) and
 * any hierarchical or local symmetries.
 * <p>
 * It also stores the parts of the structure (symmetric
 * units) involved in each axis, in addition to the way
 * to calculate them.
 * <p>
 * This is intended to provide a general axis support for
 * the multiple subunit alignment optimization and the axis
 * display in Jmol. This object is closely linked to a 
 * MultipleAlignment object that defines the symmetric units.
 * 
 * @author Aleix Lafita
 *
 */
public class SymmetryAxes {

	private List<Matrix4d> axes;
	private HashMap<Matrix4d, List<List<Integer>>> mapAxisSubunits;
	//The list of Integers defines which subunits are involved in this axis
	//and how the axis transformation is computed. The list has size 2 for 
	//the first index and any size for the second index. Describes that all
	//subunits in the first list superimposed to all subunits of the second
	//list should generate the axis transformation.

	/**
	 * Constructor. 
	 * Initializes variables only.
	 */
	public SymmetryAxes(){
		axes = new ArrayList<Matrix4d>();
		mapAxisSubunits = new HashMap<Matrix4d, List<List<Integer>>>();
	}

	/**
	 * Adds a new axis of symmetry. 
	 * The subunits that participate in this axis and their superposition
	 * relation should also be indicated.
	 * 
	 * @param axis the new axis of symmetry found
	 * @param subunits subunits participating and superposition relation
	 * 
	 * @throws IllegalArgumentException if the subunit relation is in the
	 * 			wrong format: should be double List of equal sizes.
	 */
	public void addAxis(Matrix4d axis, List<List<Integer>> subunits){

		//Check correct format of subunit relations
		if (subunits.size() != 2){
			throw new IllegalArgumentException(
					"Wrong subunit relations format: should be double List");
		} else if (subunits.get(0).size() != subunits.get(1).size()){
			throw new IllegalArgumentException(
					"Wrong subunit relations format: not equal List sizes");
		}
		axes.add(axis);
		mapAxisSubunits.put(axis, subunits);
	}

	/**
	 * Return all axes of symmetry of the structure.
	 * 
	 * @return axes of symmetry.
	 */
	public List<Matrix4d> getAxes(){
		return axes;
	}

	/**
	 * Return the subunit superposition relation needed to obtain the
	 * axis of symmetry (which subunits and in which order have to be
	 * superimposed to obtain the axis).
	 *  
	 * @param axis the symmetry axis to calculate.
	 * @return the double List of subunit relations, or null if the
	 * 			axis is not already stored.
	 */
	public List<List<Integer>> getSubunitRelation(Matrix4d axis){
		return mapAxisSubunits.get(axis);
	}


	/**
	 * Return a List of all the transformations that need to be applied
	 * to a determinate subunit in order to superimpose it to all the
	 * others. 
	 * The returned List of transformations can contain repeated matrices,
	 * because it can be that the same axis have to be applied multiple
	 * times to a single subunit.
	 * 
	 * @param subunit the subunit index
	 * @return List of transformation matrices. It can be empty
	 */
	public List<Matrix4d> getSubunitTransforms(int subunit){

		List<Matrix4d> allTransforms = new ArrayList<Matrix4d>();

		for (Matrix4d t:axes){
			List<Integer> rigid = mapAxisSubunits.get(t).get(1);
			List<Integer> moved = mapAxisSubunits.get(t).get(1);
			int m = moved.indexOf(subunit);
			int r = rigid.indexOf(subunit);
			//Handle special cases for CLOSE symmetry
			if (r==-1 || r>m){
				for (int i=0; i<m+1; i++){
					allTransforms.add(t);
				}
			}
		}
		return allTransforms;
	}

}
