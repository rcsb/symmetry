package org.biojava3.structure.dbscan;

import java.util.SortedSet;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.align.symm.CEMirrorSymm;
import org.biojava3.structure.utils.SimpleLog;
import org.rcsb.fatcat.server.PdbChainKey;

/**
 * A quick class to find mirror symmetries.
 * 
 * @author Andreas Prlic
 * @deprecated Use {@link CEMirrorSymm} instead
 */
@Deprecated
public class FindMirrorSymmetries {
	

	
	public static void main(String[] args){
		SortedSet<PdbChainKey> reps = GetRepresentatives.getRepresentatives(40);
		AtomCache cache = new AtomCache();
		
		String filename = cache.getPath()+System.getProperty("file.seperator")
				+"findMirrorSymmetries.log";
		SimpleLog.setLogFilename(filename);
		
		int total = 0;
		int symmetric = 0;
		for ( PdbChainKey r : reps){
			try {
				String name = r.toName();
				Atom[] ca1 = cache.getAtoms(name);
				Atom[] ca2 = cache.getAtoms(name);

				Atom[] ca2M = CEMirrorSymm.reverseCA2(ca2);
				CEMirrorSymm.mirrorCoordinates(ca2M);
				
				AFPChain afp = FindMirrorSymmetries.align(ca1,ca2M,name, false);
				afp.setAlgorithmName(CeMain.algorithmName);

				//boolean isSignificant =  afp.isSignificantResult();
				boolean isSignificant = false;
				if ( afp.getTMScore() > 0.3 && afp.isSignificantResult()) {
					isSignificant = true;
					//StructureAlignmentDisplay.display(afp, ca1, ca2M);
				}
				total++;
				if ( isSignificant)
					symmetric++;
				String log = name + "\t" + String.format("%.2f", afp.getProbability()) + String.format(" %.2f", afp.getTMScore());
				log += " " + symmetric + "/" + total + " (" + (symmetric/(float)total*100f)+"%)";
				
								
				SimpleLog.write(log);
				
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}
	



	/**
	 * Aligns two proteins, checking for circular permutations.
	 * 
	 * <p>When checking for mirror symmetries, create ca2 as follows:
	 * <pre>
	 * Atom[] ca2O = StructureTools.cloneCAArray(ca1);
	 * Atom[] ca2 = FindMirrorSymmetries.reverseCA2(ca2O);
	 * FindMirrorSymmetries.mirrorCoordinates(ca2);
	 * </pre>
	 * @param ca1 CA atoms of a protein
	 * @param ca2 The mirrored and reversed version of the same protein
	 * @param name the protein's name
	 * @param showMatrix If true, displays the distance matrix after masking diagonals
	 * @return
	 * @throws StructureException
	 */
	public static AFPChain align(Atom[] ca1, Atom[] ca2, String name, boolean showMatrix) throws StructureException {
		/*
		//Atom[] ca2clone = SymmetryTools.cloneAtoms(ca2);
		Atom[] ca2clone = StructureTools.duplicateCA2(ca2);
		CeParameters params = new CeParameters();
		
		CECalculator calculator = new CECalculator(params);

		//Build alignment ca1 to ca2-ca2
		AFPChain afpChain = new AFPChain();
		afpChain.setName1(name);
		afpChain.setName2(name);
		afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);

		Matrix origM = new Matrix( calculator.getMatMatrix());
		// symmetry hack, disable main diagonale

//		for ( int i = 0 ; i< rows; i++){
//			for ( int j = 0 ; j < cols ; j++){
//				int diff = Math.abs(i-j);
//
//				if ( diff < 15 ){
//					origM.set(i,j, 99);
//				}
//				int diff2 = Math.abs(i-(j-ca2.length/2));
//				if ( diff2 < 15 ){
//					origM.set(i,j, 99);
//				}
//			}
//		}
		
		
		if ( showMatrix)
			SymmetryTools.showMatrix(origM, "mirror matrix");
		
		//showMatrix(origM, "original CE matrix");
		Matrix clone =(Matrix) origM.clone();
		calculator.setMatMatrix(clone.getArray());
		//calculator.setMatMatrix(diffDistMax.getArray());
		calculator.traceFragmentMatrix( afpChain,ca1, ca2clone);
		calculator.nextStep( afpChain,ca1, ca2clone);

		//afpChain.setAlgorithmName("SpeedCE");
		//afpChain.setVersion("0.0000001");
		afpChain.setDistanceMatrix((Matrix)origM.clone());

		double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2clone);
		afpChain.setTMScore(tmScore);
				
		CeCPMain.postProcessAlignment(afpChain, ca1, ca2clone, calculator);
		
		return afpChain;
		*/
		CeMain ce = (CeMain) StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
		
		AFPChain alignment = ce.align(ca1, ca2);
		return alignment;
	}
}