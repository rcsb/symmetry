package org.biojava.nbio.structure.align.symm.gui;

import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.gui.AlignmentTextPanel;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.Block;
import org.biojava.nbio.structure.align.model.BlockImpl;
import org.biojava.nbio.structure.align.model.BlockSet;
import org.biojava.nbio.structure.align.model.BlockSetImpl;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.jama.Matrix;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JScrollPane;

/**
 * Class that provides the visualizations options for the Symmetry analysis: multiple sequence alignments,
 * multiple structural alignments, sequence panel, etc.
 * 
 * @author Aleix Lafita
 * 
 */
public class SymmetryDisplay {

	/**
	 * Method that displays two superimposed subunits as a structural alignment in jmol.
	 * Use the new method for multiple subunit alignment.
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 */
	@Deprecated
	public static void display2SuperimposedSubunits(AFPChain afpChain, Atom[] ca1, Atom[] ca2){
		
		//Create the atom arrays corresponding to the first and second subunits only
		Atom[] ca1block = Arrays.copyOfRange(ca1, afpChain.getOptAln()[0][0][0], afpChain.getOptAln()[0][0][afpChain.getOptAln()[0][0].length-1]+1);
		Atom[] ca2block = Arrays.copyOfRange(ca2, afpChain.getOptAln()[0][1][0], afpChain.getOptAln()[0][1][afpChain.getOptAln()[0][1].length-1]+1);
		
		//Modify the optimal alignment to include only one subunit (block)
		int[][][] optAln = new int[1][][];
		int[][] block = new int[2][];
		//Normalize the residues of the subunits, to be in the range of ca1block and ca2block
		int start1 = afpChain.getOptAln()[0][0][0];
		int start2 = afpChain.getOptAln()[0][1][0];
		int[] chain1 = new int[afpChain.getOptLen()[0]];
		int[] chain2 = new int[afpChain.getOptLen()[0]];
		for (int i=0; i<chain1.length; i++){
			chain1[i] = afpChain.getOptAln()[0][0][i] - start1;
			chain2[i] = afpChain.getOptAln()[0][1][i] - start2;
		}
		block[0] = chain1;
		block[1] = chain2;
		optAln[0] = block;
		
		//Modify the AFP chain to adapt the new optimal alignment of two subunits.
		AFPChain displayAFP = new AFPChain();
		
		try {
			displayAFP = AlignmentTools.replaceOptAln(optAln, displayAFP, ca1block, ca2block);
		} catch (StructureException e1) {
			e1.printStackTrace();
		}

		//Set the name of the protein
		displayAFP.setName1(afpChain.getName1()+" su1");
		displayAFP.setName2(afpChain.getName2()+" su2");
		
		//Display the alignment of the subunits
		StructureAlignmentJmol jmolPanel = new StructureAlignmentJmol(displayAFP, ca1block, ca2block);
		jmolPanel.evalString("hide ligand;");
	}
	
	/**
	 * Method that displays all superimposed subunits as a multiple alignment in jmol.
	 * @param afpChain AFP alignment with subunits segmented as blocks
	 * @param ca1
	 * @throws StructureAlignmentException 
	 * @throws StructureException 
	 * @throws IOException 
	 */
	public static void displaySuperimposedSubunits(AFPChain afpChain, Atom[] ca1) throws StructureException, StructureAlignmentException, IOException{
		
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
		
		for (int bk=0; bk<afpChain.getBlockNum(); bk++){
			
			//Normalize the residues of the subunits, to be in the range of subunit
			int start = afpChain.getOptAln()[bk][0][0];
			List<Integer> chain = new ArrayList<Integer>();
			
			for (int i=0; i<afpChain.getOptAln()[bk][0].length; i++)
				chain.add(afpChain.getOptAln()[bk][0][i] - start);
			
			block.getAlignRes().add(chain);
		}
		
		//Display the alignment of the subunits
		displayMultipleAlignment(multAln, afpChain.getBlockRotationMatrix()[0], afpChain.getBlockShiftVector()[0]);
	}
	
	/**
	 * Method that displays a multiple structure alignment of the given alignment.
	 * @param afpChain AFP alignment with subunits segmented as blocks
	 * @param rot Matrix of rotation
	 * @param shift Atom with translation coordinates
	 * @throws StructureAlignmentException 
	 * @throws StructureException 
	 * @throws IOException 
	 */
	public static void displayMultipleAlignment(MultipleAlignment multAln, Matrix rot, Atom shift) throws StructureException, StructureAlignmentException, IOException{
		
		//Apply the symmetry transformation different number of times
		int size = multAln.size();

		List<Atom[]> atomArrays = multAln.getEnsemble().getAtomArrays();
		for (int i=0; i<size; i++){
			if (atomArrays.get(i).length < 1) 
				throw new StructureException("Length of atoms arrays is too short! " + atomArrays.get(i).length);
		}
		
		List<Atom[]> rotatedAtoms = new ArrayList<Atom[]>();
		//Rotate the atom coordinates of all the structures
		for (int i=0; i<size; i++){
			Structure displayS = atomArrays.get(i)[0].getGroup().getChain().getParent().clone();
			Atom[] rotCA = StructureTools.getRepresentativeAtomArray(displayS);
			//Rotate the structure the appropiate number of times
			for (int k=0; k<i; k++){
				Calc.rotate(displayS, rot);
				Calc.shift(displayS, shift);
			}
			rotatedAtoms.add(rotCA);
		}
		MultipleAlignmentJmol jmol = new MultipleAlignmentJmol(multAln, rotatedAtoms);
		jmol.setTitle(jmol.getStructure().getPDBHeader().getTitle());
		//Include the rotation axis also in the multiple alignments
		RotationAxis axis = new RotationAxis(rot, shift);
		jmol.evalString(axis.getJmolScript(rotatedAtoms.get(0)));
	}
	
	/**
	 * Method that displays the multiple structure alignment of all symmetry rotations in jmol.
	 * @param afpChain AFP alignment with subunits segmented as blocks
	 * @param ca1
	 * @throws StructureAlignmentException 
	 * @throws StructureException 
	 * @throws IOException 
	 */
	public static void displayMultipleAlignment(AFPChain afpChain, Atom[] ca1) throws StructureException, StructureAlignmentException, IOException{
			
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
			
			for (int k=0; k<order; k++){
				List<Integer> chain = new ArrayList<Integer>();
				
				for (int i=0; i<afpChain.getOptAln()[0][0].length; i++)
					chain.add(afpChain.getOptAln()[(bk+k)%order][0][i]);
				
				block.getAlignRes().add(chain);
			}
		}
		displayMultipleAlignment(multAln, afpChain.getBlockRotationMatrix()[0], afpChain.getBlockShiftVector()[0]);
	}
	
	/**
	 * Method that displays a multiple sequence alignment of the subunits.
	 */
	public static void showAlignmentImage(AFPChain afpChain, Atom[] ca1){

		JFrame frame = new JFrame();
		
		String result = calculateMultipleAln(afpChain, ca1);

		String title = afpChain.getAlgorithmName() + " V."+afpChain.getVersion() + " : " + afpChain.getName1();
		frame.setTitle(title);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		AlignmentTextPanel txtPanel = new AlignmentTextPanel();
		txtPanel.setText(result);

		JMenuBar menu = MenuCreator.getAlignmentTextMenu(frame,txtPanel,afpChain);

		frame.setJMenuBar(menu);
		JScrollPane js = new JScrollPane();
		js.getViewport().add(txtPanel);
		js.getViewport().setBorder(null);
		//js.setViewportBorder(null);
		//js.setBorder(null);
		//js.setBackground(Color.white);

		frame.getContentPane().add(js);
		frame.pack();      
		frame.setVisible(true);
	}
	
	/**
	 * 	Method that calculates a multiple alignment of the subunits as a String to be displayed.
	 */
	private static String calculateMultipleAln(AFPChain afpChain, Atom[] ca1) {
		
		//Create the gropus of aligned residues from the afpChain. Dimensions group[order][res nr.]
		int order = afpChain.getBlockNum();
		List<List<Integer>> groups = new ArrayList<List<Integer>>();
		for (int i=0; i<order; i++){
			List<Integer> residues = new ArrayList<Integer>();
			for (int j=0; j<afpChain.getOptLen()[i]; j++){
				residues.add(afpChain.getOptAln()[i][0][j]);
			}
			groups.add(residues);
		}
		int subunitSize = groups.get(0).size();
		
		String result = new String();  //The string that stores the multiple alignment
		String[] subunits = new String[order];  //the alignment string of every subunit, with its gaps
		String[] symb = new String[order];   //the string with the alignment information (black if not aligned)
		for (int i=0; i<order; i++){
			if (i<10) subunits[i] = "Subunit 0"+(i+1)+": ";
			else subunits[i] = "Subunit "+(i+1)+": ";
		}
		Arrays.fill(symb, "            ");
		
		int[] position = new int[order];  //the positions in every subunit
		int[] next_position = new int[order];
		Arrays.fill(position, 0);
		for (int j=0; j<order; j++){
			next_position[j] = groups.get(j).get(position[j]);
		}
		
	    //Loop for every residue to see if it is included or not in the alignments and if there is a gap
		while (true){
			boolean stop = false;
			int gaps = 0;
			char[] provisional = new char[order];
			
			//If the position is higher than the subunit insert a gap
			for (int j=0; j<order; j++){
				if (position[j]>subunitSize-1){
					provisional[j] = '-';
					gaps++;
				}
				else {
					//If the next position is lower than the residue aligned there is a gap, so increment the gap
					int res = groups.get(j).get(position[j]);
					if (next_position[j]<res){
						provisional[j] = StructureTools.get1LetterCode(ca1[next_position[j]].getGroup().getPDBName());
						gaps++;
					}
					//If they are the same do not increment gap and consider a gap in case other subunits have residues in between
					else {
						provisional[j] = '-';
					}
				}
			}
			//If all sequences have gaps means that there are unaligned residues, so include them in the alignment
			if (gaps==order){
				for (int j=0; j<order; j++){
					subunits[j] += StructureTools.get1LetterCode(ca1[next_position[j]].getGroup().getPDBName());
					symb[j] += " ";
					next_position[j]++;
				}
			}
			//If there are not gaps add the aligned residues
			else if (gaps==0){
				for (int j=0; j<order; j++){
					subunits[j] += StructureTools.get1LetterCode(ca1[next_position[j]].getGroup().getPDBName());
					symb[j] += "|";
					position[j]++;
					next_position[j]++;
				}
			}
			//If only some subunits have gaps consider this information and add gaps to the subunits with missing residues
			else{
				for (int j=0; j<order; j++){
					if (provisional[j] == '-'){
						subunits[j] += '-';
					}
					else{
						subunits[j] += StructureTools.get1LetterCode(ca1[next_position[j]].getGroup().getPDBName());
						next_position[j] ++;
					}
					symb[j] += " ";
				}
			}
			//Stop if all of the subunits have been analyzed until the end (all residues in the group)
			stop = true;
			for (int q=0; q<order; q++){
				if (position[q] < subunitSize)
					stop = false;
			}
			//Stop if any subunit has reached the end of the molecule
			for (int q=0; q<order; q++){
				if (next_position[q] > ca1.length-1)
					stop = true;
			}
			if (stop) break;
	    }
		result+="\nMultiple Subunit Alignment for "+afpChain.getName1()+"\nSubunit size: "+subunitSize+"\n\n";
		for (int j=0; j<order; j++){
			result+=subunits[j]+" "+ca1[groups.get(j).get(subunitSize-1)].getGroup().getResidueNumber().getSeqNum()+":"+ca1[groups.get(j).get(subunitSize-1)].getGroup().getChainId()+"\n";
			if (j<order-1){
				result += symb[j]+"\n";
			}
		}
		
		return result;
	}
}
