package org.biojava.nbio.structure.align.symm.gui;

import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.gui.ScaleableMatrixPanel;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.utils.SymmetryTools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;
import javax.vecmath.Matrix4d;

/**
 * Class that provides the visualizations options for the 
 * Symmetry analysis: multiple sequence alignments,
 * multiple structural alignments, alignment panel, etc.
 * 
 * @author Aleix Lafita
 * 
 */
public class SymmetryDisplay {

	/**
	 * Method that displays all superimposed subunits as a multiple alignment in jmol.
	 * @param afpChain AFP alignment with subunits segmented as blocks
	 * @param ca1
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Deprecated
	public static void displaySuperimposedSubunits(AFPChain afpChain, Atom[] ca1) throws StructureException, IOException{

		//Create new structure containing the atom arrays corresponding to separate subunits
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();
		for (int i=0; i<afpChain.getBlockNum(); i++){
			Structure newStr = new StructureImpl();
			Chain newCh = new ChainImpl();
			newStr.addChain(newCh);
			Atom[] subunit = Arrays.copyOfRange(ca1, afpChain.getOptAln()[i][0][0], afpChain.getOptAln()[i][0][afpChain.getOptAln()[i][0].length-1]+1);
			for (int k=0; k<subunit.length; k++)newCh.addGroup((Group) subunit[k].getGroup().clone());
			atomArrays.add(StructureTools.getAtomCAArray(newCh));
		}

		//Initialize a new MultipleAlignment to store the aligned subunits, each one inside a BlockSet
		MultipleAlignment multAln = new MultipleAlignmentImpl();
		multAln.getEnsemble().setAtomArrays(atomArrays);
		multAln.getEnsemble().setAlgorithmName(afpChain.getAlgorithmName());
		List<String> structureNames = new ArrayList<String>();
		for(int i=0; i<afpChain.getBlockNum(); i++) structureNames.add(afpChain.getName1());
		multAln.getEnsemble().setStructureNames(structureNames);

		//All the residues are aligned in one block only
		BlockSet blockSet = new BlockSetImpl(multAln);
		Block block = new BlockImpl(blockSet);
		block.setAlignRes(new ArrayList<List<Integer>>());

		for (int bk=0; bk<afpChain.getBlockNum(); bk++){

			//Normalize the residues of the subunits, to be in the range of subunit
			int start = afpChain.getOptAln()[bk][0][0];
			List<Integer> chain = new ArrayList<Integer>();

			for (int i=0; i<afpChain.getOptAln()[bk][0].length; i++)
				chain.add(afpChain.getOptAln()[bk][0][i] - start);

			block.getAlignRes().add(chain);
		}

		//Set the transformations
		List<Matrix4d> transforms = new ArrayList<Matrix4d>();
		Matrix4d original = Calc.getTransformation(afpChain.getBlockRotationMatrix()[0], afpChain.getBlockShiftVector()[0]);
		for (int str=0; str<multAln.size(); str++){
			Matrix4d transform = (Matrix4d) original.clone();
			for (int st=0; st<str; st++) transform.mul(original);
			transforms.add(transform);
		}
		multAln.setTransformations(transforms);
		MultipleAlignmentScorer.calculateScores(multAln);

		//Display the alignment of the subunits
		MultipleAlignmentDisplay.display(multAln).evalString(new RotationAxis(afpChain).getJmolScript(ca1));

	}

	/**
	 * Method that displays the multiple structure alignment of all symmetry rotations in jmol.
	 * @param afpChain AFP alignment with subunits segmented as blocks
	 * @param ca1
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Deprecated
	public static void displayMultipleAlignment(AFPChain afpChain, Atom[] ca1) throws StructureException, IOException{

		//Create a list with multiple references to the atom array of the structure
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();
		for (int i=0; i<afpChain.getBlockNum(); i++) atomArrays.add(ca1);

		//Initialize a new MultipleAlignment to store the aligned subunits, each one inside a BlockSet
		MultipleAlignment multAln = new MultipleAlignmentImpl();
		multAln.getEnsemble().setAtomArrays(atomArrays);
		multAln.getEnsemble().setAlgorithmName(afpChain.getAlgorithmName());
		List<String> structureNames = new ArrayList<String>();
		for(int i=0; i<afpChain.getBlockNum(); i++) structureNames.add(afpChain.getName1());
		multAln.getEnsemble().setStructureNames(structureNames);
		int order = afpChain.getBlockNum();

		for (int bk=0; bk<order; bk++){

			//Every subunit has a new BlockSet
			BlockSet blockSet = new BlockSetImpl(multAln);
			Block block = new BlockImpl(blockSet);
			block.setAlignRes(new ArrayList<List<Integer>>());

			for (int k=0; k<order; k++){
				List<Integer> chain = new ArrayList<Integer>();

				for (int i=0; i<afpChain.getOptAln()[0][0].length; i++)
					chain.add(afpChain.getOptAln()[(bk+k)%order][0][i]);

				block.getAlignRes().add(chain);
			}
		}
		//Set the transformations
		List<Matrix4d> transforms = new ArrayList<Matrix4d>();
		Matrix4d original = Calc.getTransformation(afpChain.getBlockRotationMatrix()[0], afpChain.getBlockShiftVector()[0]);
		for (int str=0; str<multAln.size(); str++){
			Matrix4d transform = (Matrix4d) original.clone();
			for (int st=0; st<str; st++) transform.mul(original);
			transforms.add(transform);
		}
		multAln.setTransformations(transforms);

		//Display the multiple alignment of the rotations
		MultipleAlignmentDisplay.display(multAln).evalString(new RotationAxis(afpChain).getJmolScript(ca1));
	}

	/**
	 * Displays a multiple alignment of the symmetry subunits.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static void subunitDisplay(MultipleAlignment msa) 
			throws StructureException {

		MultipleAlignment subunits = SymmetryTools.toSubunitAlignment(msa);
		MultipleAlignmentDisplay.display(subunits);
	}

	/**
	 * Displays a multiple alignment of the whole structure transformations
	 * colored by blocks, corresponding to the subunits.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static void fullDisplay(MultipleAlignment msa) 
			throws StructureException {

		MultipleAlignment full = SymmetryTools.toFullAlignment(msa);

		MultipleAlignmentJmol jmol = MultipleAlignmentDisplay.display(full);
		jmol.setColorByBlocks(true);
	}

	/**
	 * Show a Matrix in a new JFrame.
	 * Is this method used antwhere? If so, it should use the biojava
	 * code to display matrices. - Aleix
	 * 
	 * @param m Matrix to display
	 * @param string title of the frame
	 */
	@Deprecated
	public static void showMatrix(Matrix m, String string) {
		ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
		JFrame frame = new JFrame();

		smp.setMatrix((Matrix)m.clone());
		//smp.getMatrixPanel().setScale(0.8f);

		frame.setTitle(string);
		frame.getContentPane().add(smp);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}
}
