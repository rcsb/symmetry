package org.biojava3.structure.align.symm.quarternary;

import java.util.Arrays;
import java.util.List;

public class AlignTest {

	public static void main(String[] args) {
		
		// Aspartate Transcarbamylase (ATCase)
		List<String> pdbIds = Arrays.asList(new String[]{"1Q95","1R0B","1R0C","1RAA","1RAB","1RAC"});
		
		// Transthyretin
	//	List<String> pdbIds = Arrays.asList(new String[]{"3KGT","1TTC"});
        AligQuaternaryStructure aligner = new AligQuaternaryStructure(pdbIds);
        aligner.align();
	}
}
