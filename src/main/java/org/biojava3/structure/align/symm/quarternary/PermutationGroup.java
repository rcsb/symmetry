
package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Peter
 */
public class PermutationGroup {
    List<List<Integer>> permutations = new ArrayList<List<Integer>>();

    public void addPermutation(List<Integer> permutation) {
        if (!permutations.contains(permutation)) {
            permutations.add(permutation);
        }
    }

    public List<Integer> getPermutation(int index) {
    	return permutations.get(index);
    }
    
    public int getOrder() {
    	return permutations.size();
    }
    
 
    /**
     * Ways to complete group: 
     * - combinations of permutations pi x pj
     * - combinations with itself p^k
     * 
     */
    public void completeGroup() {
    	int n = permutations.size();
        for (int i = 0; i < permutations.size(); i++) {
            for (int j = i; j < permutations.size(); j++) {
                List<Integer> p = combine(permutations.get(i), permutations.get(j));
                addPermutation(p);
            }
        }
        // repeat iteratively until no new permutation are created
        // the  following 
        if (permutations.size() > n) {
        	int m = permutations.size();
 //       	System.out.println("completeGroup iteration");
        	completeGroup();
        	if (permutations.size() > m) {
        		System.out.println("complete group iteration: " + m +"/" + permutations.size());
        	}
        }
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Permutation Group: " + permutations.size() + " permutation");
        for (List<Integer> permutation : permutations) {
            sb.append(permutation.toString());
        }
        return sb.toString();
    }

    public static List<Integer> combine(List<Integer> permutation1, List<Integer> permutation2) {
        List<Integer> intermediate = new ArrayList<Integer>(permutation1.size());
        for (int i = 0, n = permutation1.size(); i < n; i++) {
            intermediate.add(permutation2.get(permutation1.get(i)));
        }
        return intermediate;
    }

    public static int getOrder(List<Integer> permutation) {
        List<Integer> copy = new ArrayList<Integer>(permutation);
        for (int i = 0, n = permutation.size(); i < n; i++) {
            copy = combine(copy, permutation);
            if (copy.equals(permutation)) {
                return i + 1;
            }
        }
        return 0;
    }
}
