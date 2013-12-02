package org.biojava3.structure.dbscan;

import java.util.SortedSet;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.utils.SimpleLog;
import org.biojava3.structure.utils.SymmetryTools;

/**
 * 
 * @deprecated
 */
@Deprecated
public class FindRotationSymmetries {
	public static void main(String[] args){
		SortedSet<StructureName> reps = GetRepresentatives.getRepresentatives(40);
		AtomCache cache = new AtomCache();
		
		SimpleLog.setLogFilename(cache.getPath()+System.getProperty("file.seperator")
				+"findRotationSymm.log");
		
		int total = 0;
		int symmetric = 0;
		for ( StructureName r : reps){
			try {
				String name = r.getName();
				Atom[] ca1 = cache.getAtoms(name);
				Atom[] ca2 = cache.getAtoms(name);

				AFPChain afp = FindRotationSymmetries.align(ca1,ca2,name,false);
				afp.setAlgorithmName(CeMain.algorithmName);

				//boolean isSignificant =  afp.isSignificantResult();
				boolean isSignificant = false;
				if ( afp.getTMScore() > 0.3 && afp.isSignificantResult())
					isSignificant = true;
				total++;
				if ( isSignificant)
					symmetric++;
				String log = (name + "\t" + String.format("%.2f", afp.getProbability()) + String.format(" %.2f", afp.getTMScore()));
				log += (" " + symmetric + "/" + total + " (" + (symmetric/(float)total*100f)+"%)");
				
				SimpleLog.write(log);
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}


	/**
	 * Detects off-diagonal similarity between two proteins. Typically this is
	 * used with ca1==ca2 to find rotational symmetry in the protein.
	 * @param ca1 The first protein
	 * @param ca2 The second protein, typically the same as ca1
	 * @param name A name to be used for the aligment
	 * @param showMatrix If true, displays the distance matrix after masking diagonals
	 * @return
	 * @throws StructureException
	 */
	public static AFPChain align(Atom[] ca1, Atom[] ca2, String name, boolean showMatrix) throws StructureException {
		// duplicate ca2
		Atom[] ca2m = StructureTools.duplicateCA2(ca2);
		int rows = ca1.length ;
		int cols = ca2m.length ; // duplicated length

		CeParameters params = new CeParameters();
		
		//params.setMaxNrIterationsForOptimization(0);
		CECalculator calculator = new CECalculator(params);

		//Build alignment ca1 to ca2-ca2
		AFPChain afpChain = new AFPChain();
		afpChain.setName1(name);
		afpChain.setName2(name);
		afpChain = calculator.extractFragments(afpChain, ca1, ca2m);

		Matrix origM = new Matrix( calculator.getMatMatrix());
		// symmetry hack, disable main diagonale

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				if ( diff < 15 ){
					origM.set(i,j, 99);
				}
				int diff2 = Math.abs(i-(j-ca2m.length/2));
				if ( diff2 < 15 ){
					origM.set(i,j, 99);
				}
			}
		}

		
		if ( showMatrix)
			SymmetryTools.showMatrix(origM, "Rotating Matrix, main diagonal disabled");
		
		
		//showMatrix(origM, "original CE matrix");
		calculator.setMatMatrix(((Matrix)origM.clone()).getArray());
		//calculator.setMatMatrix(diffDistMax.getArray());
		calculator.traceFragmentMatrix( afpChain,ca1, ca2m);
		calculator.nextStep( afpChain,ca1, ca2m);

		afpChain.setAlgorithmName("CE-rotation " + name);
		afpChain.setVersion("0.0000001");
		
		// reset the distance matrix, since this is modified by nextStep
		afpChain.setDistanceMatrix(origM);

		double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2m);
		afpChain.setTMScore(tmScore);

		afpChain = CeCPMain.postProcessAlignment(afpChain, ca1, ca2m, calculator);
		
		return afpChain;
	}
}
