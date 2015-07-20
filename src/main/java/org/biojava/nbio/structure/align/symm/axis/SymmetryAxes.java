package org.biojava.nbio.structure.align.symm.axis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Matrix4d;

/**
 * Data Structure that stores all the symmetry axis that describe 
 * the symmetry of a structure. Generalizes to all types of symmetry, 
 * the classic ones (Cn, Dn) and any hierarchical or local symmetries.
 * <p>
 * It also stores the parts of the structure (symmetric units) involved 
 * in each axis, in addition to the way to calculate them.
 * <p>
 * This is intended to provide a general axis support for the multiple 
 * subunit alignment optimization and the axis display in Jmol. This 
 * object is related to a MultipleAlignment object that defines the 
 * symmetric units.
 * 
 * @author Aleix Lafita
 *
 */
public class SymmetryAxes {

	/**
	 * List of all symmetry axis. They are sorted from higher to lower 
	 * in the symmetry hierarchy, where higher means that they apply 
	 * more globally and lower means that they apply to a local region 
	 * of the higher axis division.
	 */
	private List<Matrix4d> axes;

	/**
	 * Tree of the symmetry axes hierarchy of a structure. Every index
	 * is the axis index and the integer stored means the number of parts
	 * this axis divides the structure in.
	 */
	private List<Integer> tree;

	/**
	 * Matrix of size [subunits][axes]. The first index of the matrix 
	 * indicates which subunit is considered, and the second index
	 * indicates which axis. The integer stored means how many times
	 * the transformation has to be applied to the subunit, and a value
	 * of 0 means that this subunit is not affected by the axis.
	 */
	private List<List<Integer>> subunitTransforms;

	/**
	 * The list of integers defines which subunits are involved in this axis
	 * and how the axis transformation is computed. The list has size 2 for 
	 * the first index and any size for the second index. Describes that all
	 * subunits in the first list superimposed to all subunits of the second
	 * list should generate the axis transformation.
	 */
	private Map<Integer, List<List<Integer>>> mapAxisSubunits;

	/**
	 * Constructor. 
	 * Initializes variables only.
	 */
	public SymmetryAxes(){
		axes = new ArrayList<Matrix4d>();
		mapAxisSubunits = new HashMap<Integer, List<List<Integer>>>();
		subunitTransforms = new ArrayList<List<Integer>>();
		tree = new ArrayList<Integer>();
	}

	/**
	 * Adds a new axis of symmetry. 
	 * The subunits that participate in this axis and their superposition
	 * relation should also be indicated.
	 * 
	 * @param axis the new axis of symmetry found
	 * @param superposition subunits participating and superposition relation
	 * @param subunits number of times the transformation is applied to every
	 * 			subunit. index1=subunit, index2=times.
	 * @param division number of parts that this axis divides the structure in
	 * 
	 * @throws IllegalArgumentException if the subunit relation is in the
	 * 			wrong format: should be double List of equal sizes.
	 */
	public void addAxis(Matrix4d axis, List<List<Integer>> superposition, 
			List<Integer> subunits, Integer division) {

		//Check correct format of subunit relations
		if (superposition.size() != 2){
			throw new IllegalArgumentException(
					"Wrong superposition format: should be double List.");
		} else if (superposition.get(0).size() != superposition.get(1).size()){
			throw new IllegalArgumentException(
					"Wrong superposition format: not equal List sizes.");
		}
		if (division < 2) throw new IllegalArgumentException(
				"A symmetry axis should divide a structure in > 2 parts");

		axes.add(axis);
		tree.add(division);

		//Extend the double List by the necessary rows
		while (subunitTransforms.size() < subunits.size()){
			List<Integer> list = new ArrayList<Integer>();
			for (int a=0; a<axes.size()-1; a++){
				list.add(0);
			}
			subunitTransforms.add(list);
		}
		for (int su=0; su<subunitTransforms.size(); su++){
			Integer nTimes = 0;
			if (su<subunits.size()){
				nTimes = subunits.get(su);
			}
			subunitTransforms.get(su).add(nTimes);
		}

		//Check that the subunit indices match
		for (int c=0; c<2; c++){
			for (int p=0; p<superposition.get(c).size(); p++){
				Integer su = superposition.get(c).get(p);
				if (su>subunitTransforms.size()){
					throw new IllegalArgumentException(
							"Subunit indicex in superposition out of bounds");
				}
			}
		}
		mapAxisSubunits.put(axes.size()-1, superposition);
	}

	/**
	 * Updates an axis of symmetry, after the superposition changed.
	 * 
	 * @param index old axis index
	 * @param newAxis
	 */
	public void updateAxis(Integer index, Matrix4d newAxis){
		axes.set(index, newAxis);
	}

	/**
	 * Return all elementary axes of symmetry of the structure, that is,
	 * the axes stored in the List as unique and from which all the symmetry
	 * axes are constructed.
	 * 
	 * @return axes elementary axes of symmetry.
	 */
	public List<Matrix4d> getElementaryAxes(){
		return axes;
	}

	/**
	 * Return the subunit superposition relation needed to obtain the
	 * axis of symmetry (which subunits and in which order have to be
	 * superimposed to obtain the axis).
	 *  
	 * @param index the axis index
	 * @return the double List of subunit relations, or null if the
	 * 			axis is not already stored.
	 */
	public List<List<Integer>> getSubunitRelation(Integer index){
		return mapAxisSubunits.get(index);
	}


	/**
	 * Return the transformation that needs to be applied to a determinate 
	 * subunit so that all get superimposed to the same point. 
	 * 
	 * @param subunit the subunit index
	 * @return transformation matrix for the subunit
	 */
	public Matrix4d getSubunitTransform(int subunit){

		List<Matrix4d> allTransforms = new ArrayList<Matrix4d>();
		Matrix4d transform = new Matrix4d();
		transform.setIdentity();

		for (int a=0; a<axes.size(); a++){
			Matrix4d t = (Matrix4d) axes.get(a).clone();
			Matrix4d clone = (Matrix4d) t.clone();
			//Pack the Matrices when they are equal
			for (int i=1; i<subunitTransforms.get(subunit).get(a); i++){
				clone.mul(t);
			}
			if (subunitTransforms.get(subunit).get(a)>0) allTransforms.add(clone);
		}
		//Multiply the matrices in the inverse order as they have to be applied
		//for (int t=0; t<allTransforms.size(); t++){
		for (int t=allTransforms.size()-1; t>=0; t--){
			transform.mul(allTransforms.get(t));
		}

		return transform;
	}

	/**
	 * Return all symmetry axes of of the structure, the set of axes that
	 * describe all parts of the structure. This combines the elementary
	 * axes using the symmetry hierarchy tree and returns a bigger set.
	 * Use this method to display the axes.
	 * 
	 * @return axes all symmetry axes of the structure.
	 */
	public List<Matrix4d> getSymmetryAxes(){

		/*List<Matrix4d> symmAxes = new ArrayList<Matrix4d>();
		TODO
		//For every elementary axis do
		for (int t=1; t<axes.size(); t++){
			symmAxes.add((Matrix4d) axes.get(t).clone());
			
			//Consider all its parent axes combinations
			for (int i=0; i<t; i++){

				//Means how many parents to consider
				for (int n=i; n<t; n++){
					Matrix4d axis = (Matrix4d) axes.get(t).clone();
					Matrix4d parent = axes.get(n);
					
					for (int k=0; k<tree.get(i); k++){
						Matrix4d clone = new Matrix4d();
						clone.setIdentity();

						for (int j=0; j<k; j++){
							clone.mul(parent);
						}
						Matrix4d invert = (Matrix4d) clone.clone();
						invert.invert();
						clone.mul(axis);
						clone.mul(invert);
						axis = clone;
					}
				}
			}
		}
		return symmAxes;*/
		return getElementaryAxes();
	}
}
