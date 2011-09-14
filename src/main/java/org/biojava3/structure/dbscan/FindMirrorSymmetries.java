package org.biojava3.structure.dbscan;

import java.util.SortedSet;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;

import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.utils.SimpleLog;
import org.biojava3.structure.utils.SymmetryTools;
import org.rcsb.fatcat.server.PdbChainKey;

//TODO Should be @Deprecated and made a CE-Symm subclass
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

				Atom[] ca2M = reverseCA2(ca2);
				FindMirrorSymmetries.mirrorCoordinates(ca2M);
				
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
	 * Reverses an array of atoms.
	 * Really only useful for the detection of mirror symmetries, of which there
	 * are only few in the PDB.
	 * 
	 * @param ca2 Array to be reversed
	 * @return A cloned and reversed copy of ca2
	 * @throws StructureException
	 */
	public static Atom[] reverseCA2(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length];

		int pos = ca2clone.length - 1;

		Chain c = new ChainImpl();
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos--;
		}


//		// Duplicate ca2!
//		for (Atom a : ca2){
//			Group g = (Group)a.getGroup().clone();
//			c.addGroup(g);
//			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);
//
//			pos--;
//		}

		return ca2clone;


	}
	
	/**
	 * Creates a mirror image of a structure along the X axis.
	 * 
	 * @param ca2O The array of atoms to be modified
	 */
	public static void mirrorCoordinates(Atom[] ca2O) {
		for(int i=0;i<ca2O.length;i++) {
			//ca2O[i].setX(-ca2O[i].getX());
			Group g = ca2O[i].getGroup();
			for ( Atom a : g.getAtoms()){
				a.setX(-a.getX());
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
	}
}