package org.biojava.nbio.structure.align.symm.gui;

import java.awt.event.KeyEvent;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.MyAlignmentLoadListener;
import org.biojava.nbio.structure.align.gui.MyDistMaxListener;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 *  Create a menu for the Symmetry analysis GUI. Adapted from MenuCreator in biojava.
 *	
 *	@author lafita
 *
 */
public class SymmetryMenu extends MenuCreator {
	
	//Menu Options for the Symmetry Display
	public static final String SUBUNIT_DISPLAY = "Subunit Superimposition";
	public static final String SUBUNIT_ALIGN = "Multiple Subunit Alignment";
	public static final String SYMMETRY = "New Symmetry Analysis";
	
	/** 
	 *  Provide a JMenuBar that can be added to a JFrame
	 */
	public static JMenuBar initMenu(JFrame frame, SymmetryJmol parent, AFPChain afpChain){

		JMenuBar menu = new JMenuBar();

		//FILE tab
		JMenu file= new JMenu("File");
		file.setMnemonic(KeyEvent.VK_F);
		file.getAccessibleContext().setAccessibleDescription("File Menu");

		if ( parent != null){
			JMenuItem loadF = getLoadMenuItem();
			loadF.addActionListener(new MyAlignmentLoadListener(parent));
			file.add(loadF);
		}

		JMenuItem saveF = getSaveAlignmentMenuItem(afpChain);
		file.add(saveF);

		JMenuItem openPDB = getShowPDBMenuItem();
		file.add(openPDB);

		JMenuItem openI = getOpenPDBMenuItem();

		file.add(openI);

		if ( parent != null){
			JMenuItem exportI =  getExportPDBMenuItem(parent);

			file.add(exportI);
		}

		JMenuItem openDBI = getDBResultMenuItem();
		file.add(openDBI);
		file.addSeparator();

		if ( parent != null){
			JMenuItem print = getPrintMenuItem();
			print.addActionListener(parent.getJmolPanel());

			file.add(print);
		}
		file.addSeparator();

		JMenuItem closeI = getCloseMenuItem(frame);

		file.add(closeI);

		JMenuItem exitI = getExitMenuItem();
		file.add(exitI);
		menu.add(file);

		//ALIGN tab
		JMenu align = new JMenu("Align");
		align.setMnemonic(KeyEvent.VK_A);
		JMenuItem pairI = getPairwiseAlignmentMenuItem();
		align.add(pairI);

		menu.add(align);

		//VIEW tab
		JMenu view = new JMenu("View");
		view.getAccessibleContext().setAccessibleDescription("View Menu");
		view.setMnemonic(KeyEvent.VK_V);

		if ( parent != null){
			JMenuItem aligpI = MenuCreator.getIcon(parent,ALIGNMENT_PANEL);
			aligpI.setMnemonic(KeyEvent.VK_L);
			aligpI.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_M, keyMask));
			view.add(aligpI);

			JMenuItem textI = MenuCreator.getIcon(parent,TEXT_ONLY);
			textI.setMnemonic(KeyEvent.VK_T);
			view.add(textI);

			JMenuItem pairsI = MenuCreator.getIcon(parent,PAIRS_ONLY);
			pairsI.setMnemonic(KeyEvent.VK_P);
			view.add(pairsI);

			JMenuItem textF = MenuCreator.getIcon(parent,FATCAT_TEXT);
			textF.setMnemonic(KeyEvent.VK_F);
			view.add(textF);
		}

		if ( afpChain != null){
			JMenuItem distMax = new  JMenuItem("Show Distance Matrices");
			distMax.setMnemonic(KeyEvent.VK_D);
			distMax.addActionListener(new MyDistMaxListener(afpChain));
			view.add(distMax);

			JMenuItem dotplot = new JMenuItem("Show Dot Plot");
			dotplot.setMnemonic(KeyEvent.VK_O);
			dotplot.addActionListener(new DotPlotListener(afpChain));
			view.add(dotplot);
		}

		menu.add(view);
		
		//SYMMETRY tab
		JMenu sym = new JMenu("Symmetry");
		sym.setMnemonic(KeyEvent.VK_S);
		
		JMenuItem subunits = new JMenuItem(SUBUNIT_DISPLAY);
		subunits.addActionListener(parent);
		subunits.setMnemonic(KeyEvent.VK_D);
		
		JMenuItem mulAln = new JMenuItem(SUBUNIT_ALIGN);
		mulAln.addActionListener(parent);
		mulAln.setMnemonic(KeyEvent.VK_M);
		
		JMenuItem newSym = new JMenuItem(SYMMETRY);
		newSym.addActionListener(parent);
		newSym.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Y, keyMask));
		newSym.setMnemonic(KeyEvent.VK_Y);
		
		sym.add(subunits);
		sym.add(mulAln);
		sym.add(newSym);

		menu.add(sym);

		//HELP tab
		JMenu about = new JMenu("Help");
		about.setMnemonic(KeyEvent.VK_H);

		JMenuItem helpM = getHelpMenuItem();
		about.add(helpM);

		JMenuItem aboutM = getAboutMenuItem();
		about.add(aboutM);

		menu.add(Box.createGlue());
		menu.add(about);

		return menu;

	}
	
	/** 
	 *  Provide a display for a symmetry analysis
	 */
	protected static void showSymmDialog(){
		SymmetryGui gui =  SymmetryGui.getInstance();
		gui.setVisible(true);
	}
}