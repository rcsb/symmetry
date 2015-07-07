package org.biojava.nbio.structure.align.symm.gui;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.MultipleStructureAligner;
import org.biojava.nbio.structure.align.gui.AlignmentCalculationRunnable;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** 
 * Extension from the Biojava class that works with one 
 * structure only, using the Symmetry specific GUI.
 *  
 * @author Aleix Lafita
 */

public class SymmetryCalc implements AlignmentCalculationRunnable {
	
	private static final Logger logger = 
			LoggerFactory.getLogger(SymmetryCalc.class);
	
	private boolean interrupted = false;
	
	private String name;
	private Structure structure;
	private SymmetryGui parent;

	/** Requests for a structure to analyze.
	 */
	public SymmetryCalc(SymmetryGui p, Structure s, String n) {
		parent = p;
		structure = s;
		name = n;
	}

	@Override
	public void run() {

		//The structure has been downloaded, now calculate the alignment ...
		MultipleStructureAligner algorithm = parent.getStructureAlignment();
		CESymmParameters params = (CESymmParameters) algorithm.getParameters();
		
		try {

			List<Atom[]> atoms = new ArrayList<Atom[]>();
			atoms.add(StructureTools.getRepresentativeAtomArray(structure));
			MultipleAlignment msa = algorithm.align(atoms);

			List<String> names = new ArrayList<String>();
			for (int su=0; su<msa.size(); su++){
				names.add(name+"_"+(su+1));
			}
			msa.getEnsemble().setStructureNames(names);

			SymmetryJmol jmol = new SymmetryJmol(msa);
			String title = jmol.getTitle();
			
			if (params != null) 
				title += " | OrderDetector=" + params.getOrderDetectorMethod()+
				" Refiner: "+params.getRefineMethod();
			jmol.setTitle(title);

		} catch (StructureException e){
			e.printStackTrace();
			logger.warn(e.getMessage());
		}
		parent.notifyCalcFinished();
	}
	
	@Override
	public void interrupt() {
		interrupted = true;
	}

	@Override
	public void cleanup() {

		parent.notifyCalcFinished();
		parent = null;
		structure = null;
	}
	
	@Override
	public void setNrCPUs(int useNrCPUs) {
		// TODO Auto-generated method stub
	}
}
