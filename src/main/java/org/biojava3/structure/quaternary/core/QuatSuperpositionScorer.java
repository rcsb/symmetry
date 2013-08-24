
package org.biojava3.structure.quaternary.core;

import java.util.List;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

/**
 *
 * @author Peter
 */
public class QuatSuperpositionScorer {

    /**
     * Returns minimum, mean, and maximum RMSD and TM-Score for two superimposed sets of subunits
     * 
     * TM score: Yang Zhang and Jeffrey Skolnick, PROTEINS: Structure, Function, and Bioinformatics 57:702–710 (2004)
     * @param subunits subunits to be scored
     * @param transformation transformation matrix
     * @param permutations permutation that determines which subunits are superposed
     * @return
     */
    public static QuatSymmetryScores calcScores(Subunits subunits, Matrix4d transformation, List<Integer> permutation) {
    	QuatSymmetryScores scores = new QuatSymmetryScores();
    
    	double minTm = Double.MAX_VALUE;
    	double maxTm = Double.MIN_VALUE;
    	double minRmsd = Double.MAX_VALUE;
    	double maxRmsd = Double.MIN_VALUE;
    	
        double totalSumTm = 0;
        double totalSumDsq = 0;
    	double totalLength = 0;
        
        Point3d t = new Point3d();
        List<Point3d[]> traces = subunits.getTraces();
        
        // loop over the Calpha atoms of all subunits
        for (int i = 0; i < traces.size(); i++) {
        	// in helical systems not all permutations involve all subunit. -1 indicates subunits that should not be permuted.
        	 if (permutation.get(i) == -1) {
             	continue;
             }
        	// get original subunit
            Point3d[] orig = traces.get(i);
            totalLength += orig.length;
            
            // get permuted subunit
            Point3d[] perm = traces.get(permutation.get(i));
            
            // calculate TM specific parameters
            int tmLen = Math.max(orig.length, 17);  // don't let d0 get negative with short sequences
            double d0 = 1.24 * Math.cbrt(tmLen - 15.0) - 1.8;
            double d0Sq = d0 * d0;

            double sumTm = 0;
            double sumDsq = 0;
            for (int j = 0; j < orig.length; j++) {
            	// transform coordinates of the permuted subunit
                t.set(perm[j]);
                transformation.transform(t);
                
                double dSq = orig[j].distanceSquared(t);
                sumTm += 1.0/(1.0 + dSq/d0Sq);
                sumDsq += dSq;
            }
            
            // scores for individual subunits
            double sTm = sumTm/tmLen;
            minTm = Math.min(minTm, sTm);
            maxTm = Math.max(maxTm, sTm);
            
            double sRmsd = Math.sqrt(sumDsq/orig.length);
            minRmsd = Math.min(minRmsd,  sRmsd);
            maxRmsd = Math.max(maxRmsd,  sRmsd);
            
            totalSumTm += sumTm;
            totalSumDsq += sumDsq;
        }
        
        // save scores for individual subunits
        scores.setMinRmsd(minRmsd);
        scores.setMaxRmsd(maxRmsd);
        scores.setMinTm(minTm);
        scores.setMaxTm(maxTm);
        
        // save mean scores over all subunits
        scores.setTm(totalSumTm/totalLength);
        scores.setRmsd(Math.sqrt(totalSumDsq/totalLength));
        
		return scores;
    }
    
}
