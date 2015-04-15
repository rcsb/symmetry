package org.biojava.nbio.structure.align.symm.gui;

import java.awt.Color;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.gui.AlignmentCalculationRunnable;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.jcolorbrewer.ColorBrewer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** Extension from the Biojava class that allows a constructor with one structure only.
 *  
 * @author lafita
 */

public class SymmetryCalc implements AlignmentCalculationRunnable {
	
	private static final Logger logger = LoggerFactory.getLogger(SymmetryCalc.class);

	boolean interrupted = false;

	String pdb;
	String name;
	Structure structure1;
	Structure structure2;
	SymmetryGui parent;

	/** Requests for a structure to analyze.
	 */
	public SymmetryCalc(SymmetryGui p, Structure s1, Structure s2, String n) {
		parent = p;
		structure1 = s1;
		structure2 = s2;
		name = n;
	}

	@Override
	public void run() {

		// the structure has been downloaded, now calculate the alignment ...

		StructureAlignment algorithm = parent.getStructureAlignment();  
		//StructurePairAligner aligner = new StructurePairAligner();
		//aligner.setDebug(true);
		try {

			Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
			Atom[] ca2 = StructureTools.getAtomCAArray(structure2);

			//System.out.println("ca1 size:" + ca1.length + " ca2 size: " + ca2.length);
			AFPChain afpChain = algorithm.align(ca1, ca2);

			afpChain.setName1(name);
			afpChain.setName2(name);
			
			CESymmParameters params = (CESymmParameters) algorithm.getParameters();
			//Set the color of the subunits
			Color[] subunitColors = null;
			CESymmParameters.SubunitColors COLOR = params.getSubunitColors();
			switch(COLOR){
			case COLOR_SET: 
				subunitColors = ColorBrewer.Set1.getColorPalette(afpChain.getBlockNum());
				break;
			case SPECTRAL:
				subunitColors = ColorBrewer.Spectral.getColorPalette(afpChain.getBlockNum());
				break;
			case PAIRED:
				subunitColors = ColorBrewer.Paired.getColorPalette(afpChain.getBlockNum());
			case GRADUAL:
				break;
			}

			SymmetryJmol jmol = SymmetryDisplay.display(afpChain, ca1, ca2, subunitColors);

			String title = jmol.getTitle();
			
			/*if ( params != null) //The title is very large for CeSymm, so disabled.
				title += " " + algorithm.getParameters().toString();*/
			jmol.setTitle(title);

			//DisplaySymmAFP.showAlignmentImage(afpChain,ca1,ca2,jmol);

			System.out.println(afpChain.toCE(ca1,ca2));

		} catch (StructureException e){
			e.printStackTrace();
			logger.warn(e.getMessage());

		}
		//logger.info("done!");

		parent.notifyCalcFinished();

	}

	/** stops what is currently happening and does not continue
	 * 
	 *
	 */
	@Override
	public void interrupt() {
		interrupted = true;
	}

	@Override
	public void cleanup() {

		parent.notifyCalcFinished();

		parent=null;
		// cleanup...

		structure1 = null;
		structure2 = null;

	}

	/** does not do anything here...
	 * 
	 */
	@Override
	public void setNrCPUs(int useNrCPUs) {
		// TODO Auto-generated method stub
		// 
	}
}