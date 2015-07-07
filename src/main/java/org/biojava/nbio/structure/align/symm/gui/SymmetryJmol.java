package org.biojava.nbio.structure.align.symm.gui;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JTextField;
import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.JmolPanel;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.gui.jmol.MyJmolStatusListener;
import org.biojava.nbio.structure.align.gui.jmol.RasmolCommandListener;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.align.webstart.AligUIManager;
import org.biojava.nbio.structure.jama.Matrix;
import org.jcolorbrewer.ColorBrewer;

/** 
 * A class that provides a simple GUI for symmetry alignments in Jmol.
 * It only displays the structure with the colored subunits and rotation axis.
 * Provides links to other Jmol displays through the menus.
 * 
 * @author Aleix Lafita
 * 
 */
public class SymmetryJmol extends AbstractAlignmentJmol {

	private Color[] subunitColors;
	private MultipleAlignment msa;
	private Atom[] ca;

	/**
	 * Empty Constructor.
	 * @throws StructureException
	 */
	public SymmetryJmol() throws StructureException {
		this(null, null);
	}

	/**
	 * Constructor with a MultipleAlignment of the subunits. Default axis.
	 * The atoms in the have to be the complete structure for all subunits,
	 * but the alignment should only contain one Block with every
	 * subunit as a new row.
	 * 
	 * @param MultipleAlignment subunit multiple alignment 
	 * (generated from optimization)
	 * @throws StructureException
	 */
	public SymmetryJmol(MultipleAlignment alignment) throws StructureException{

		this(alignment, new ArrayList<RotationAxis>());

		//Get the DEFAULT rotation axis if the user does not input it
		if (alignment.getTransformations() != null){
			Matrix4d transform = alignment.getTransformations().get(1);
			RotationAxis axis = new RotationAxis(transform);
			evalString(axis.getJmolScript(ca));
		}
	}

	/**
	 * Main Constructor with a MultipleAlignment of the subunits.
	 * The atoms in the alignment must be the complete structure,
	 * but the alignment should only contain one Block with every
	 * subunit as a new row.
	 * 
	 * @param MultipleAlignment subunit multiple alignment 
	 * (generated from optimization)
	 * @param axis set of rotation axis that describe the symmetry of the structure
	 * @throws StructureException
	 */
	public SymmetryJmol(MultipleAlignment alignment, List<RotationAxis> axis) throws StructureException {

		AligUIManager.setLookAndFeel();

		nrOpenWindows++;
		jmolPanel = new JmolPanel();
		frame = new JFrame();
		JMenuBar menu = SymmetryMenu.initJmolMenu(frame, this, alignment);
		frame.setJMenuBar(menu);

		this.msa = alignment;
		this.ca = msa.getEnsemble().getAtomArrays().get(0);
		this.subunitColors = ColorBrewer.Spectral.getColorPalette(alignment.size());

		frame.addWindowListener(new WindowAdapter()
		{
			@Override
			public void windowClosing(WindowEvent e) {

				nrOpenWindows--;
				destroy();
				if ( nrOpenWindows > 0){
					frame.dispose();
				}
				else  {
					// check if AlignmentGUI is visible..
					SymmetryGui gui = SymmetryGui.getInstanceNoVisibilityChange();
					if ( gui.isVisible()) {
						frame.dispose();
						gui.requestFocus();
					} else {
						System.exit(0);
					}
				}
			}
		});

		Container contentPane = frame.getContentPane();
		Box vBox = Box.createVerticalBox();

		jmolPanel.addMouseMotionListener(this);
		jmolPanel.addMouseListener(this);

		jmolPanel.setPreferredSize(new Dimension(DEFAULT_WIDTH,DEFAULT_HEIGHT));
		vBox.add(jmolPanel);

		// USER SCRIPTING COMMAND
		JTextField field = new JTextField();

		field.setMaximumSize(new Dimension(Short.MAX_VALUE,30));   
		field.setText(COMMAND_LINE_HELP);
		RasmolCommandListener listener = new RasmolCommandListener(jmolPanel,field) ;

		field.addActionListener(listener);
		field.addMouseListener(listener);
		field.addKeyListener(listener);
		vBox.add(field);


		/// COMBO BOXES 
		Box hBox1 = Box.createHorizontalBox();
		hBox1.add(Box.createGlue());

		String[] styles = new String[] { "Cartoon", "Backbone", "CPK", "Ball and Stick", "Ligands","Ligands and Pocket"};
		JComboBox style = new JComboBox(styles);
		hBox1.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		hBox1.add(new JLabel("Style"));
		hBox1.add(style);
		vBox.add(hBox1);
		contentPane.add(vBox);

		style.addActionListener(jmolPanel);

		String[] colorModes = new String[] { "Secondary Structure", "By Chain", "Rainbow", "By Element", "By Amino Acid", "Hydrophobicity" ,"Suggest Domains" , "Show SCOP Domains"};
		JComboBox colors = new JComboBox(colorModes);
		colors.addActionListener(jmolPanel);
		hBox1.add(Box.createGlue());
		hBox1.add(new JLabel("Color"));
		hBox1.add(colors);

		String[] colorPattelete = new String[] {"Color Set", "Spectral", "2Colors", "3Colors", "Pastel", "Reds", "Blues" ,"Greens"};
		JComboBox pattelete = new JComboBox(colorPattelete);

		pattelete.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				JComboBox source = (JComboBox) e.getSource();
				String value = source.getSelectedItem().toString();
				evalString("save selection; select *; color grey; select ligand; color CPK;");
				if (value=="Color Set"){
					subunitColors = ColorBrewer.Set1.getColorPalette(msa.size());
					colorPalette = ColorBrewer.Set1;
				} else if (value=="Spectral"){
					subunitColors = ColorBrewer.Spectral.getColorPalette(msa.size());
					colorPalette = ColorBrewer.Spectral;
				} else if (value=="2Colors"){
					subunitColors = ColorBrewer.Set1.getColorPalette(2);
					colorPalette = ColorBrewer.Set1;
				} else if (value=="3Colors"){
					subunitColors = ColorBrewer.Set1.getColorPalette(3);
					colorPalette = ColorBrewer.Set1;
				} else if (value=="Pastel"){
					subunitColors = ColorBrewer.Pastel1.getColorPalette(msa.size());
					colorPalette = ColorBrewer.Pastel1;
				} else if (value=="Reds"){
					subunitColors = ColorBrewer.Reds.getColorPalette(msa.size());
					colorPalette = ColorBrewer.Reds;
				} else if (value=="Blues"){
					subunitColors = ColorBrewer.Blues.getColorPalette(msa.size());
					colorPalette = ColorBrewer.Blues;
				} else if (value=="Greens"){
					subunitColors = ColorBrewer.Greens.getColorPalette(msa.size());
					colorPalette = ColorBrewer.Greens;
				}
				evalString(getJmolString(msa, ca, subunitColors)+"; restore selection;");
			}
		});

		hBox1.add(Box.createGlue());
		hBox1.add(new JLabel("Symmetry"));
		hBox1.add(pattelete);


		// CHeck boxes
		Box hBox2 = Box.createHorizontalBox();
		hBox2.setMaximumSize(new Dimension(Short.MAX_VALUE,30));

		JButton resetDisplay = new JButton("Reset Display");

		resetDisplay.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				System.out.println("reset!!");
				jmolPanel.executeCmd("restore STATE state_1");

			}
		});

		hBox2.add(resetDisplay); 
		hBox2.add(Box.createGlue());

		JCheckBox toggleSelection = new JCheckBox("Show Selection");
		toggleSelection.addItemListener(
				new ItemListener() {
					@Override
					public void itemStateChanged(ItemEvent e) {
						boolean showSelection = (e.getStateChange() == ItemEvent.SELECTED);

						if (showSelection){
							jmolPanel.executeCmd("set display selected");
						} else {
							jmolPanel.executeCmd("set display off");
						}
					}
				}
				);

		hBox2.add(toggleSelection);
		hBox2.add(Box.createGlue());
		vBox.add(hBox2);	


		// STATUS DISPLAY

		Box hBox = Box.createHorizontalBox();

		status = new JTextField();		
		status.setBackground(Color.white);
		status.setEditable(false);
		status.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		status.setPreferredSize(new Dimension(DEFAULT_WIDTH / 2,30));
		status.setMinimumSize(new Dimension(DEFAULT_WIDTH / 2,30));
		hBox.add(status);      
		text = new JTextField();
		text.setBackground(Color.white);
		text.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		text.setPreferredSize(new Dimension(DEFAULT_WIDTH / 2,30));
		text.setMinimumSize(new Dimension(DEFAULT_WIDTH / 2,30));
		text.setText("Display of Atom info");
		text.setEditable(false);
		hBox.add(text);
		vBox.add(hBox);

		contentPane.add(vBox);
		MyJmolStatusListener li = (MyJmolStatusListener) jmolPanel.getStatusListener();
		li.setTextField(status);
		frame.pack();
		frame.setVisible(true);

		// init coordinates
		initCoords();
		printSymmetryAxis(axis);
		resetDisplay();
	}

	private void printSymmetryAxis(List<RotationAxis> symmetryAxis){

		for (int a=0; a<symmetryAxis.size(); a++){
			String script = symmetryAxis.get(a).getJmolScript(ca, a);
			evalString(script);
		}
	}

	/**
	 * Override the action listeners of the menu items to add the new options
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		String cmd = e.getActionCommand();		    
		if (cmd.equals(MenuCreator.ALIGNMENT_PANEL)){
			if (msa == null) {
				System.err.println("Currently not displaying a symmetry!");
				return;
			}
			try {
				MultipleAlignmentDisplay.showMultipleAligmentPanel(msa, this, colorPalette);
			} catch (StructureException e1) {
				e1.printStackTrace();
			}

		} else if (cmd.equals(MenuCreator.FATCAT_TEXT)){
			if (msa == null) {
				System.err.println("Currently not displaying a symmetry!");
				return;
			}
			String result = MultipleAlignmentWriter.toFatCat(msa);
			result += AFPChain.newline;
			result += MultipleAlignmentWriter.toTransformMatrices(msa);
			MultipleAlignmentDisplay.showAlignmentImage(msa, result);

		} else if (cmd.equals(SymmetryMenu.SUBUNIT_DISPLAY)){
			if (msa == null) {
				System.err.println("Currently not displaying a symmetry!");
				return;
			}
			try {
				SymmetryDisplay.subunitDisplay(msa);
			} catch (Exception e1) {
				e1.printStackTrace();
			}

		} else if (cmd.equals(SymmetryMenu.MULTIPLE_STRUCT)){
			if (msa == null) {
				System.err.println("Currently not displaying a symmetry!");
				return;
			}
			try {
				SymmetryDisplay.fullDisplay(msa);
			} catch (Exception e1){
				e1.printStackTrace();
			}

		} else if (cmd.equals(SymmetryMenu.SYMMETRY)){

			SymmetryMenu.showSymmDialog();
		}
	}

	public static String getJmolString(MultipleAlignment msa, Atom[] ca, Color[] subunitColors) {

		StringWriter jmol = new StringWriter();
		jmol.append(DEFAULT_SCRIPT);

		List<List<Integer>> alignRes = msa.getBlocks().get(0).getAlignRes();

		for(int str=0; str < alignRes.size(); str++) {

			printJmolScript4Block(ca, alignRes, jmol, str, subunitColors);
			jmol.append("backbone 0.6 ;");
		}

		jmol.append(LIGAND_DISPLAY_SCRIPT);
		//System.out.println(jmol);
		return jmol.toString();

	}

	private static void printJmolScript4Block(Atom[] ca, List<List<Integer>> alignRes, StringWriter jmol, int str, Color[] colors) {

		Color c1 = colors[str%colors.length];

		List<String> pdb = new ArrayList<String>();
		for (int i=0; i< alignRes.get(str).size(); i++) {
			Integer pos = alignRes.get(str).get(i);
			if (pos != null) pdb.add(JmolTools.getPdbInfo(ca[pos]));
		}

		// and now select the aligned residues...
		StringBuffer buf = new StringBuffer("select ");
		int count = 0;
		for (String res : pdb){
			if (count > 0) buf.append(",");
			buf.append(res);
			count++;
		}
		buf.append("; color [" + c1.getRed() +"," + c1.getGreen() +"," +c1.getBlue()+"];");

		jmol.append(buf);
	}

	@Override
	protected void initCoords() {
		if (ca == null ){
			if (structure != null) setStructure(structure);
		} else {
			structure = ca[0].getGroup().getChain().getParent();
		}
		jmolPanel.setStructure(structure);
	}

	@Override
	public void resetDisplay() {

		if (msa != null && ca != null) {
			String script = getJmolString(msa, ca, subunitColors);
			//System.out.println(script);
			script += "select ligand; color CPK;";
			evalString(script);
			jmolPanel.evalString("save STATE state_1");
		}
	}

	@Override
	public List<Matrix> getDistanceMatrices() {
		if (msa==null) return null;
		else return Arrays.asList(msa.getEnsemble().getDistanceMatrix().get(0));
	}
}
