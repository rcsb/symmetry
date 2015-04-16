package org.biojava.nbio.structure.align.symm.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JScrollPane;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.AlignmentTextPanel;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.aligpanel.StatusDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.gui.aligpanel.SymmAligPanel;
import org.biojava.nbio.structure.align.util.AlignmentTools;

/** A utility class for visualistion of symmetry alignments (rotation of the second structure).
 * 	Adapted from the DisplayAFP class of biojava.
 * 
 * @author lafita
 *
 */
public class DisplaySymmAFP extends DisplayAFP {
	
	/** 
	 * Rotation of ca2, hetatoms2 and nucleotides2 will be done here.
	 */
	public static SymmetryJmol display(AFPChain afpChain,Group[] twistedGroups, Atom[] ca1, Atom[] ca2, List<Group> hetatms, List<Group> hetatms2, Color[] subunitColors) throws StructureException{

		List<Atom> twistedAs = new ArrayList<Atom>();

		for ( Group g: twistedGroups){
			if ( g == null )
				continue;
			if ( g.size() < 1)
				continue;
			Atom a = g.getAtom(0);
			twistedAs.add(a);
		}
		Atom[] twistedAtoms = (Atom[])twistedAs.toArray(new Atom[twistedAs.size()]);

		Atom[] arr1 = getAtomArray(ca1, hetatms);
		Atom[] arr2 = getAtomArray(twistedAtoms, hetatms2);

		// 

		//if ( hetatms2.size() > 0)
			//	System.out.println("atom after:" + hetatms2.get(0).getAtom(0));

		//if ( hetatms2.size() > 0)
		//	System.out.println("atom after:" + hetatms2.get(0).getAtom(0));

		String title =  afpChain.getAlgorithmName() + " V." +afpChain.getVersion() + " : " + afpChain.getName1();

		//System.out.println(artificial.toPDB());

		SymmetryJmol jmol = new SymmetryJmol(afpChain,arr1,arr2, subunitColors);
		jmol.setTitle(title);
		
		return jmol;
	}
	
	
	/**
	 * Method that displays two superimposed subunits in jmol.
	 * 
	 * INPUT: an AFP alignment and the protein.
	 * OUTPUT: a JmolPanel with only one subunit superimposed.
	 */
	public static void displaySuperimposedSubunits(AFPChain afpChain, Atom[] ca1, Atom[] ca2){
		
		//Create the atom arrays corresponding to the first and second subunits only
		Atom[] ca1block = new Atom[afpChain.getOptLen()[0]];
		Atom[] ca2block = new Atom[afpChain.getOptLen()[0]];
		ca1block = Arrays.copyOfRange(ca1, 0, afpChain.getOptAln()[0][0][afpChain.getOptAln()[0][0].length-1]+1);
		ca2block = Arrays.copyOfRange(ca2, afpChain.getOptAln()[0][1][0], afpChain.getOptAln()[0][1][afpChain.getOptAln()[0][1].length-1]+1);
		
		//Modify the optimal alignment to include only one subunit (block)
		int[][][] optAln = new int[1][2][afpChain.getOptLen()[0]];
		int[][] block = afpChain.getOptAln()[0];
		//Normalize the residues of the second subunit, to be in the range of ca2block
		int start = block[1][0];
		for (int i=0; i<block[1].length; i++){
			block[1][i] -= start;
		}
		optAln[0] = block;
		int[] optLens = new int[1];
		optLens[0]=optAln[0][0].length;
		
		//Modify the AFP chain to adapt the new optimal alignment of two subunits.
		AFPChain displayAFP = new AFPChain();
		try {
			displayAFP = AlignmentTools.replaceOptAln(optAln, afpChain, ca1block, ca2block);
		} catch (StructureException e1) {
			e1.printStackTrace();
		}
		
		//Another array to display is created only with the residues of the second subunit, because all (first and second, are needed to superimpose, but only the second is relevant in the alignment)
		//DOES NOT WORK, because the second subunit is not colored
		//Atom[] ca2blockDisplay = Arrays.copyOfRange(ca2, afpChain.getOptAln()[0][1][0], afpChain.getOptAln()[0][1][afpChain.getOptAln()[0][1].length-1]+1);

		//Set the name of the protein
		displayAFP.setName1(afpChain.getName1()+" su1");
		displayAFP.setName2(afpChain.getName2()+" su2");
		displayAFP.setAlgorithmName("jFatCat_rigid");  //Set the name to FatCat, because it is not a CE-Symm alignment.
		
		try {
			//Display the AFP alignment of the subunits
			StructureAlignmentJmol jmolPanel;
			jmolPanel = StructureAlignmentDisplay.display(displayAFP, ca1block, ca2block);
			jmolPanel.evalString("hide ligand;");
			
		} catch (StructureException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Method that displays a multiple alignment of the subunits.
	 */
	public static void showMulAlnImage(AFPChain afpChain, Atom[] ca1, Atom[] ca2){

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
	
	public static void showAlignmentImage(AFPChain afpChain, Atom[] ca1, Atom[] ca2, SymmetryJmol jmol, Color[] subunitColors) {
		String result = afpChain.toFatcat(ca1, ca2);

		//String rot = afpChain.toRotMat();
		//DisplayAFP.showAlignmentImage(afpChain, result + AFPChain.newline + rot);

		SymmAligPanel me = new SymmAligPanel();
		me.setSubunitColors(subunitColors);
		me.setStructureAlignmentJmol(jmol);
		me.setAFPChain(afpChain);

		JFrame frame = new JFrame();

		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);		
		frame.setTitle(afpChain.getName1() + " vs. " + afpChain.getName2() + " | " + afpChain.getAlgorithmName() + " V. " + afpChain.getVersion());
		me.setPreferredSize(new Dimension(me.getCoordManager().getPreferredWidth() , me.getCoordManager().getPreferredHeight()));

		JMenuBar menu = MenuCreator.getAlignmentTextMenu(frame,me,afpChain);
		frame.setJMenuBar(menu);

		JScrollPane scroll = new JScrollPane(me);
		scroll.setAutoscrolls(true);

		StatusDisplay status = new StatusDisplay();
		status.setAfpChain(afpChain);

		status.setCa1(ca1);
		status.setCa2(ca2);
		me.setCa1(ca1);
		me.setCa2(ca2);
		me.addAlignmentPositionListener(status);


		Box vBox = Box.createVerticalBox();
		vBox.add(scroll);
		vBox.add(status);


		frame.getContentPane().add(vBox);

		frame.pack();
		frame.setVisible(true);
		// make sure they get cleaned up correctly:
			frame.addWindowListener(me);
			frame.addWindowListener(status);
	}
	
	/**
	 * 	Method that calculates a multiple alignment of the subunits as a String to be displayed.
	 */
	public static String calculateMultipleAln(AFPChain afpChain, Atom[] ca1) {
		
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
			subunits[i] = "Subunit "+(i+1)+": ";
		}
		Arrays.fill(symb, "           ");
		
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
			result+=subunits[j]+"\n";
			if (j<order-1){
				result += symb[j]+"\n";
			}
		}
		
		return result;
	}
}
