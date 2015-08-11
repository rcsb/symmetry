package org.biojava.nbio.structure.align.symm.benchmark.external;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;

/**
 * Finds parts of alignments that use the main diagonal.
 * @author dmyerstu
 */
public class TrivialIdentityFinder {

	public static void main(String[] args) {
		if (args.length != 1) {
			System.err.println("Usage: " + TrivialIdentityFinder.class.getSimpleName() + " directory-of-fasta-files");
			return;
		}
		File dir = new File(args[0]);
		print(dir.listFiles());
	}

	public static void print(File... fastaFiles) {
		Map<File,Integer> map = new HashMap<File,Integer>();
		for (File file : fastaFiles) {
			int count = -1;
			try {
				count = countResiduesFromFasta(file);
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
			map.put(file, count);
		}
		System.out.println("=========================================");
		System.out.println("Found " + fastaFiles.length + " files");
		System.out.println("=========================================");
		for (Map.Entry<File,Integer> entry : map.entrySet()) {
			System.out.println(entry.getKey().getName() + "\t" + entry.getValue());
		}
		System.out.println("=========================================");
	}
	
	public static int countResiduesFromFasta(File fastaFile) throws StructureException, IOException {
		return countResidues(SymDFasta.getAlignment(fastaFile));
	}
	
	public static int countResidues(AFPChain afpChain) throws StructureException {
		return countResidues(AlignmentTools.alignmentAsMap(afpChain));
	}

	public static int countResidues(Map<Integer, Integer> alignmentAsMap) {
		int count = 0;
		for (Map.Entry<Integer,Integer> entry : alignmentAsMap.entrySet()) {
			if (entry.getKey().equals(entry.getValue())) count++;
		}
		return count;
	}

}
