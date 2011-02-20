package org.biojava3.structure.dbscan;

import java.util.SortedSet;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;
import org.rcsb.fatcat.server.PdbChainKey;

public class FindRotationSymmetries {

	public static void main(String[] args){
		SortedSet<PdbChainKey> reps = GetRepresentatives.getRepresentatives();
		AtomCache cache = new AtomCache("/Users/andreas/WORK/PDB/",true);
		FindRotationSymmetries me = new FindRotationSymmetries();

		int total = 0;
		int symmetric = 0;
		for ( PdbChainKey r : reps){
			try {
				String name = r.toName();
				Atom[] ca1 = cache.getAtoms(name);
				Atom[] ca2 = cache.getAtoms(name);

				AFPChain afp = me.align(ca1,ca2,name);
				afp.setAlgorithmName(CeMain.algorithmName);

				//boolean isSignificant =  afp.isSignificantResult();
				boolean isSignificant = false;
				if ( afp.getTMScore() > 0.3 && afp.isSignificantResult())
					isSignificant = true;
				total++;
				if ( isSignificant)
					symmetric++;
				System.out.print(name + "\t" + String.format("%.2f", afp.getProbability()) + String.format(" %.2f", afp.getTMScore()));
				System.out.println(" " + symmetric + "/" + total + " (" + (symmetric/(float)total*100f)+"%)");
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}


	private AFPChain align(Atom[] ca1, Atom[] ca2, String name) throws StructureException {
		int rows = ca1.length ;
		int cols = ca2.length ;

		CeParameters params = new CeParameters();
		params.setCheckCircular(true);
		CECalculator calculator = new CECalculator(params);

		//Build alignment ca1 to ca2-ca2
		AFPChain afpChain = new AFPChain();
		afpChain.setName1(name);
		afpChain.setName2(name);
		afpChain = calculator.extractFragments(afpChain, ca1, ca2);

		Matrix origM = new Matrix( calculator.getMatMatrix());
		// symmetry hack, disable main diagonale

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				if ( diff < 15 ){
					origM.set(i,j, 99);
				}
				int diff2 = Math.abs(i-(j-ca2.length/2));
				if ( diff2 < 15 ){
					origM.set(i,j, 99);
				}
			}
		}

		//showMatrix(origM, "original CE matrix");
		calculator.setMatMatrix(origM.getArray());
		//calculator.setMatMatrix(diffDistMax.getArray());
		calculator.traceFragmentMatrix( afpChain,ca1, ca2);
		calculator.nextStep( afpChain,ca1, ca2);

		afpChain.setAlgorithmName("CE-rotation");
		afpChain.setVersion("0.0000001");


		double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
		afpChain.setTMScore(tmScore);

		return afpChain;


	}
}
