package org.biojava.nbio.structure.align.symm.census2.utils;

import java.io.IOException;
import java.util.Map;

import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.symm.census2.Results;
import org.biojava.nbio.structure.align.symm.census2.Significance;
import org.biojava.nbio.structure.align.symm.census2.SignificanceFactory;

/**
 * Find results in a census that use the trivial alignment.
 * @author dmyersturnbull
 */
public class TrivialAlignmentFinder {

	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + CensusDifferences.class.getSimpleName() + " census-1.xml");
			return;
		}
		find(Results.fromXML(args[0]));
	}

	public static void find(Results census) {
		Significance sig = SignificanceFactory.forCeSymmOrd();
		System.out.println("----");
		int count = 0;
		for (Result r : census.getData()) {
			Map<Integer, Integer> map = r.getAlignmentMapping().getSimpleFunction();
			int nBad = 0;
			StringBuilder bad = new StringBuilder();
			for (Map.Entry<Integer, Integer> entry : map.entrySet()) {
				if (entry.getKey() == entry.getValue()) {
					nBad++;
					bad.append(entry.getKey() + " ");
				}
			}
			if (nBad > 3) {
				count++;
				System.out.println(r.getScopId() + "\t" + bad);
				if (sig.isSignificant(r)) {
					System.out.println("*****SIGNIFICANT*****");
					System.out.println();
				}
//				System.out.println(AlignmentTools.toConciseAlignmentString(r.getAlignmentMapping().getSimpleFunction()));
//				System.out.println();
			}
		}
		System.out.println(count);
	}

}
