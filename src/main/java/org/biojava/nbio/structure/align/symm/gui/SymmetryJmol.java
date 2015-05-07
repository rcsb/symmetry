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
import java.util.List;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JTextField;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.JmolPanel;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.gui.jmol.MyJmolStatusListener;
import org.biojava.nbio.structure.align.gui.jmol.RasmolCommandListener;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.align.webstart.AligUIManager;
import org.biojava.nbio.structure.gui.util.color.ColorUtils;

/** 
 * A class that provides a simple GUI for symmetry alignments in Jmol.
 * It only displays the structure with the colored subunits and rotation axis.
 * Provides links to other Jmol displays (Muliple Alignments) through the menus.
 * 
 * @author lafita
 * 
 */
public class SymmetryJmol extends AbstractAlignmentJmol {
	   
	Color[] subunitColors;
	AFPChain afpChain;
	Atom[] ca;
	
	/**
	 * Empty Constructor.
	 * @throws StructureException
	 */
	public SymmetryJmol() throws StructureException{
		this(null,null,null);
	}
	
	/**
	 * Constructor without color specifications.
	 * @param afpChain
	 * @param ca1
	 * @throws StructureException
	 */
	public SymmetryJmol(AFPChain afpChain, Atom[] ca1) throws StructureException {
		this(afpChain,ca1,null);
	}
	
	/**
	 * Main Constructor with all the arguments.
	 * @param afpChain
	 * @param ca1
	 * @param subunitColors
	 * @throws StructureException
	 */
	public SymmetryJmol(AFPChain afpChain, Atom[] ca1, Color[] subunitColors) throws StructureException {
		  
	      AligUIManager.setLookAndFeel();

	      nrOpenWindows++;
	      jmolPanel = new JmolPanel();
	      frame = new JFrame();
	      JMenuBar menu = SymmetryMenu.initMenu(frame,this, afpChain);
	      frame.setJMenuBar(menu);
	      
	      this.afpChain = afpChain;
	      this.ca = ca1;
	      this.subunitColors = subunitColors;

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
	      resetDisplay();
	      
		  //String with the options to format the Symmetry Starting view
		  RotationAxis axis = new RotationAxis(afpChain);
		  String cmd = axis.getJmolScript(ca1);
		  cmd += "select ligand; color grey;";
		  jmolPanel.evalString(cmd);
		
	}
	
	/**
	 * Override the action listeners of the menu items to add the new options
	 */
	@Override
   	public void actionPerformed(ActionEvent e) {
		String cmd = e.getActionCommand();		    
		if (cmd.equals(SymmetryMenu.SEQUENCE_PANEL)){
		    if ( afpChain == null) {
		       System.err.println("Currently not viewing an alignment!");
		       return;
		    }
		    SymmetryDisplay.showSequenceImage(afpChain, ca, this, subunitColors);
		    
		} else if (cmd.equals(SymmetryMenu.SUBUNIT_DISPLAY)){
	    	 if ( afpChain == null) {
	              System.err.println("Currently not viewing a symmetry!");
	              return;
	    	 }
	         try {
					SymmetryDisplay.displaySuperimposedSubunits(afpChain, ca);
			} catch (Exception e1) {
				e1.printStackTrace();
			}
	         
		} else if (cmd.equals(SymmetryMenu.SEQUENCE_ALIGN)){
	    	 if (afpChain == null) {
	              System.err.println("Currently not viewing a symmetry!");
	              return;
	    	 }
	         try {
					SymmetryDisplay.showAlignmentImage(afpChain, ca);
			} catch (Exception e1) {
				e1.printStackTrace();
			}
	         
		} else if (cmd.equals(SymmetryMenu.MULTIPLE_STRUCT)){
	    	 if ( afpChain == null) {
	              System.err.println("Currently not viewing a symmetry!");
	              return;
	          }
	    	  try {
				SymmetryDisplay.displayMultipleAlignment(afpChain, ca);
			} catch (Exception e1){
				e1.printStackTrace();
			}
	    	  
		} else if (cmd.equals(SymmetryMenu.SYMMETRY)){
			
	    	  SymmetryMenu.showSymmDialog();
	    }
      }
   
	   public static String getJmolScript4Block(AFPChain afpChain, Atom[] ca, int blockNr, Color[] subunitColors){
		   int blockNum = afpChain.getBlockNum();
		   
		   if ( blockNr >= blockNum)
			   return DEFAULT_SCRIPT;
			   		   
		   int[] optLen = afpChain.getOptLen();
		   int[][][] optAln = afpChain.getOptAln();
	
		   if ( optLen == null)
			   return DEFAULT_SCRIPT;
	
		   StringWriter jmol = new StringWriter();
		   jmol.append(DEFAULT_SCRIPT);
		      
		   printJmolScript4Block(ca, blockNum, optLen, optAln, jmol, blockNr, subunitColors);
		   
		   jmol.append(LIGAND_DISPLAY_SCRIPT);
		   //System.out.println(jmol);
		   return jmol.toString();
	
		   
	   }
	   
	
	   private static String getJmolString(AFPChain afpChain, Atom[] ca1, Color[] subunitColors)
	   {
	
	      int blockNum = afpChain.getBlockNum();      
	      int[] optLen = afpChain.getOptLen();
	      int[][][] optAln = afpChain.getOptAln();
	
	      if ( optLen == null)
	         return DEFAULT_SCRIPT;
	
	      StringWriter jmol = new StringWriter();
	      jmol.append(DEFAULT_SCRIPT);
	      
	      for(int bk = 0; bk < blockNum; bk ++)       {
	
	         printJmolScript4Block(ca1, blockNum, optLen, optAln, jmol, bk, subunitColors);
	      }
	      
	      jmol.append(LIGAND_DISPLAY_SCRIPT);
	      //System.out.println(jmol);
	      return jmol.toString();
	      
	   }
	   
	   private static void printJmolScript4Block(Atom[] ca, int blockNum,
				int[] optLen, int[][][] optAln, StringWriter jmol, int bk, Color[] colors) {
			//the block nr determines the color...
			 int colorPos = bk;
			 
			 Color c1;
			 //If the colors for the block are specified in AFPChain use them, otherwise the default ones are calculated
			 if (colors==null){
				 
				 if ( colorPos > ColorUtils.colorWheel.length){
				    colorPos = ColorUtils.colorWheel.length % colorPos ;
				 }
				 
				 Color end1 = ColorUtils.rotateHue(ColorUtils.cyan,    (1.0f  / 24.0f) * (blockNum +1)  ) ;
				 	 
				 c1   = ColorUtils.getIntermediate(ColorUtils.cyan, end1, blockNum, bk);
			 }
			 else {
				 int n = colors.length;
				 c1   = colors[colorPos%n];
			 }
			 
			 List<String> pdb1 = new ArrayList<String>();
			 for ( int i=0;i< optLen[bk];i++) {
			    int pos1 = optAln[bk][0][i];
			    pdb1.add(JmolTools.getPdbInfo(ca[pos1]));
			 }

			 // and now select the aligned residues...
			 StringBuffer buf = new StringBuffer("select ");
			 int count = 0;
			 for (String res : pdb1 ){
			    if ( count > 0)
			       buf.append(",");
			    buf.append(res);
			    count++;
			 }
			 buf.append("; backbone 0.6 ; color [" + c1.getRed() +"," + c1.getGreen() +"," +c1.getBlue()+"];");
			 
			 //buf.append("; set display selected;");
			 // now color this block:
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
		
		if (afpChain != null && ca != null) {
	         String script = getJmolString(afpChain,ca,subunitColors);
	         System.out.println(script);
	         evalString(script);
	         jmolPanel.evalString("save STATE state_1");
	      }
		
	}
}
