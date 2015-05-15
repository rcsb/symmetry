package org.biojava.nbio.structure.align.symm.gui;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.gui.AlignmentCalculationRunnable;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** Extension from the Biojava class that allows a constructor with one structure only.
 *  
 * @author Aleix Lafita
 */

public class SymmetryCalc implements AlignmentCalculationRunnable {
	
	private static final Logger logger = LoggerFactory.getLogger(SymmetryCalc.class);

	boolean interrupted = false;

	String pdb;
	String name;
	Structure structure;
	SymmetryGui parent;

	/** Requests for a structure to analyze.
	 */
	public SymmetryCalc(SymmetryGui p, Structure s, String n) {
		parent = p;
		structure = s;
		name = n;
	}

	@Override
	public void run() {

		// the structure has been downloaded, now calculate the alignment ...
		StructureAlignment algorithm = parent.getStructureAlignment();
		CESymmParameters params = (CESymmParameters) algorithm.getParameters();
		
		try {

			Atom[] ca1 = StructureTools.getAtomCAArray(structure);
			Atom[] ca2 = StructureTools.getAtomCAArray(structure.clone());

			//System.out.println("ca1 size:" + ca1.length + " ca2 size: " + ca2.length);
			AFPChain afpChain = algorithm.align(ca1, ca2);

			afpChain.setName1(name);
			afpChain.setName2(name);

			SymmetryJmol jmol = new SymmetryJmol(afpChain, ca1);

			String title = jmol.getTitle();
			
			if (params != null) title += " | OrderDetector=" + params.getOrderDetectorMethod()+" Refiner: "+params.getRefineMethod();
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

		structure = null;
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