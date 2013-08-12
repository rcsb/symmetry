
package org.biojava3.structure.quaternary.core;

import java.util.List;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

/**
 *
 * @author Peter
 */
public class QuatSuperpositionScorer {
    private Subunits subunits = null;

    public QuatSuperpositionScorer(Subunits subunits) {
        this.subunits = subunits;
    }
    
    public double calcCalphaRMSD(Matrix4d transformation, List<Integer> permutation) {
    	// rmsd can only be calculated if the permuted subunits
    	// are identical in sequence
    	if (! hasEquivalentSubunits(permutation)) {
    		return -1.0;
    	}
    	
        int len = 0;
        double distanceSq = 0;
        
        // calculate total, tertiary, and quat. rmsd.
        
        Point3d t = new Point3d();
        List<Point3d[]> traces = subunits.getTraces();

        for (int i = 0; i < traces.size(); i++) {
            Point3d[] orig = traces.get(i);
            len += orig.length;
            // transform each CA trace
            Point3d[] perm = traces.get(permutation.get(i));
            
            for (int j = 0; j < perm.length; j++) {
                t.set(perm[j]);
                transformation.transform(t);
                distanceSq += orig[j].distanceSquared(t);
            }
        }
        
        return Math.sqrt(distanceSq/len);
    }
    
    /**
     * Returns the TM-Score for two superimposed sets of coordinates
     * Yang Zhang and Jeffrey Skolnick, PROTEINS: Structure, Function, and Bioinformatics 57:702–710 (2004)
     * @param transformation transformation matrix
     * @param permutations permutation that determines which subunits are superposed
     * @param lengthNative total length of native sequence
     * @return
     */
    public double calcCalphaMinTMScore(Matrix4d transformation, List<Integer> permutation) {
    	// rmsd can only be calculated if the permuted subunits
    	// are identical in sequence
    	if (! hasEquivalentSubunits(permutation)) {
    		return -1.0;
    	}
    	
    	double tmScoreMin = Double.MAX_VALUE;
        
        Point3d t = new Point3d();
        List<Point3d[]> traces = subunits.getTraces();
        
        for (int i = 0; i < traces.size(); i++) {
            Point3d[] orig = traces.get(i);
            
            // transform each CA trace
            Point3d[] perm = traces.get(permutation.get(i));
            
            // what should be the length of native and aligned sequence here??
            double d0 = 1.24 * Math.cbrt(orig.length - 15.0) - 1.8;
            double d0Sq = d0 * d0;

            double tmScore = 0;
            for (int j = 0; j < orig.length; j++) {
                t.set(perm[j]);
                transformation.transform(t);
                tmScore += 1.0/(1.0 + orig[j].distanceSquared(t)/d0Sq);
            }
            tmScore /= orig.length;
            tmScoreMin = Math.min(tmScoreMin,  tmScore);
        }
        
        return tmScoreMin;
    }
    
    /**
     * Returns the TM-Score for two superimposed sets of coordinates
     * Yang Zhang and Jeffrey Skolnick, PROTEINS: Structure, Function, and Bioinformatics 57:702–710 (2004)
     * @param transformation transformation matrix
     * @param permutations permutation that determines which subunits are superposed
     * @param lengthNative total length of native sequence
     * @return
     */
    public double calcCalphaWtTMScore(Matrix4d transformation, List<Integer> permutation) {
    	// rmsd can only be calculated if the permuted subunits
    	// are identical in sequence
    	if (! hasEquivalentSubunits(permutation)) {
    		return -1.0;
    	}
    	
    	double tmScoreWt = 0;
    	double totalLength = 0;
        
        Point3d t = new Point3d();
        List<Point3d[]> traces = subunits.getTraces();
        for (int i = 0; i < traces.size(); i++) {
            totalLength += traces.get(i).length;
        }

        for (int i = 0; i < traces.size(); i++) {
            Point3d[] orig = traces.get(i);
            
            // transform each CA trace
            Point3d[] perm = traces.get(permutation.get(i));
            
            // what should be the length of native and aligned sequence here??
            double d0 = 1.24 * Math.cbrt(orig.length - 15.0) - 1.8;
            double d0Sq = d0 * d0;

            for (int j = 0; j < orig.length; j++) {
                t.set(perm[j]);
                transformation.transform(t);
                tmScoreWt += 1.0/(1.0 + orig[j].distanceSquared(t)/d0Sq);
            }
        }
        
        return tmScoreWt/totalLength;
    }
    
    /**
     * Returns true if the specified permutation permutes 
     * equivalent subunits. 
     * 
     * @param permutation
     * @return
     */
    private boolean hasEquivalentSubunits(List<Integer> permutation) {
    	List<Integer> sequenceClusterIds = subunits.getSequenceClusterIds();
  
    	for (int i = 0; i < sequenceClusterIds.size(); i++) {
    		int j = permutation.get(i);
    		if (sequenceClusterIds.get(i) != sequenceClusterIds.get(j)) {
    			return false;
    		}
    	}
    	return true;
    }

}
