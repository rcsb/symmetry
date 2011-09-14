/**
 * 
 */
package org.biojava3.structure.align.symm;

import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;

/**
 * @author Spencer Bliven
 *
 */
public class CEMirrorSymm extends CeMain {
	static String algorithmName = "jCE-mirror-symmetry";
	static String version = "0.1";
	
	/**
	 * Aligns the first protein against the second in reverse topological order.
	 * This allows the detection of mirror symmetries.
	 * 
	 * @param ca1 The first protein
	 * @param ca2 The second protein, typically a clone of the first protein
	 * @param param A {@link CeParameters} object
	 * @see org.biojava3.structure.align.symm.CeSymm#align(org.biojava.bio.structure.Atom[], org.biojava.bio.structure.Atom[], java.lang.Object)
	 */
	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param)
			throws StructureException {
		
		// Optionally, mirror the coordinates
		// TODO preference for this
		mirrorCoordinates(ca2);

		Atom[] ca2m = reverseCA2(ca2);

		AFPChain afpChain = super.align(ca1, ca2m, param);
		
		// try to set name2, which was lost in the clone
		try {
			afpChain.setName2(ca2[0].getGroup().getChain().getParent().getName());
		} catch( Exception e) {}
		
		postProcessAlignment(afpChain);
		
		return afpChain;
	}

	private void postProcessAlignment(AFPChain afpChain) {
		// reverse the optAln
		reverseOptAln(afpChain);
		
		// reverse the distance matrix
		afpChain.setDistanceMatrix(reverseMatrixCols(afpChain.getDistanceMatrix()));
		
		// reverse the ca2 matrix
		Matrix distMat2 = afpChain.getDisTable2();
		distMat2 = reverseMatrixRows(distMat2);
		distMat2 = reverseMatrixCols(distMat2);
		afpChain.setDisTable2(distMat2);
		
	}
	
	private static void reverseOptAln(AFPChain afpChain) {
		int ca2len = afpChain.getCa2Length();
		
		int[][][] optAln = afpChain.getOptAln();
		int[] optLen = afpChain.getOptLen();
		
		for(int block = 0; block<afpChain.getBlockNum(); block++) {
			for(int pos = 0; pos< optLen[block]; pos++) {
				optAln[block][1][pos] = ca2len -1 - optAln[block][1][pos]; 
			}
		}
		
		afpChain.setOptAln(optAln);
	}
	
	private static Matrix reverseMatrixRows(Matrix mat) {
		int[] reversed = new int[mat.getRowDimension()];
		for(int i=0;i<reversed.length;i++) {
			reversed[i] = reversed.length-i-1;
		}
		
		Matrix revMat = mat.getMatrix(reversed, 0, mat.getColumnDimension()-1);
		return revMat;
	}

	private static Matrix reverseMatrixCols(Matrix mat) {
		int[] reversed = new int[mat.getColumnDimension()];
		for(int i=0;i<reversed.length;i++) {
			reversed[i] = reversed.length-i-1;
		}
		
		Matrix revMat = mat.getMatrix(0, mat.getRowDimension()-1, reversed);
		return revMat;
	}
	
	/**
	 * @see org.biojava3.structure.align.symm.CeSymm#getAlgorithmName()
	 */
	@Override
	public String getAlgorithmName() {
		return CEMirrorSymm.algorithmName;
	}

	/**
	 * @see org.biojava3.structure.align.symm.CeSymm#getVersion()
	 */
	@Override
	public String getVersion() {
		return CEMirrorSymm.version; 
	}

	/**
	 * Reverses an array of atoms.
	 * Really only useful for the detection of mirror symmetries, of which there
	 * are only few in the PDB.
	 * 
	 * @param ca2 Array to be reversed
	 * @return A cloned and reversed copy of ca2
	 * @throws StructureException
	 */
	public static Atom[] reverseCA2(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length];

		int pos = ca2clone.length - 1;

		Chain c = new ChainImpl();
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos--;
		}


//		// Duplicate ca2!
//		for (Atom a : ca2){
//			Group g = (Group)a.getGroup().clone();
//			c.addGroup(g);
//			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);
//
//			pos--;
//		}

		return ca2clone;


	}
	
	/**
	 * Creates a mirror image of a structure along the X axis.
	 * 
	 * @param ca2O The array of atoms to be modified
	 */
	public static void mirrorCoordinates(Atom[] ca2O) {
		for(int i=0;i<ca2O.length;i++) {
			//ca2O[i].setX(-ca2O[i].getX());
			Group g = ca2O[i].getGroup();
			for ( Atom a : g.getAtoms()){
				a.setX(-a.getX());
			}
		}
	}
	
	public static void main(String[] args) {
		String name = "1qys";
		
		try {
			StructureAlignmentFactory.addAlgorithm(new CEMirrorSymm());
			CEMirrorSymm ce = (CEMirrorSymm) StructureAlignmentFactory.getAlgorithm(CEMirrorSymm.algorithmName);
			
			AtomCache cache = new AtomCache();
			
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
						
			AFPChain afpChain = ce.align(ca1, ca2);
			
			StructureAlignmentDisplay.display(afpChain, ca1, ca2);

		} catch (StructureException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
