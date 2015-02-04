package org.biojava.nbio.structure.align.symm.census2.utils;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.symm.census2.Results;
import org.biojava.nbio.structure.align.symm.census2.Significance;
import org.biojava.nbio.structure.align.symm.census2.SignificanceFactory;
import org.biojava.nbio.structure.align.symm.census3.stats.CensusStatUtils;

/**
 * Find differences between two census files.
 * @author dmyersturnbull
 */
public class CensusDifferences {

	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: " + CensusDifferences.class.getSimpleName() + " census-1.xml census-2.xml");
			return;
		}
		printDifferences(Results.fromXML(args[0]), Results.fromXML(args[1]));
	}

	private static void printDifferences(Results census1, Results census2) {

		Significance sig = SignificanceFactory.forCeSymmOrd();
		
		Map<String, Float> m1 = new HashMap<String, Float>();
		for (Result r : census1.getData()) {
			m1.put(r.getScopId(), r.getAlignment().getTmScore());
		}
		Map<String, Integer> a1 = new HashMap<String, Integer>();
		for (Result r : census1.getData()) {
			a1.put(r.getScopId(), r.getAlignment().getAlignLength());
		}

		Map<String, Float> m2 = new HashMap<String, Float>();
		for (Result r : census2.getData()) {
			m2.put(r.getScopId(), r.getAlignment().getTmScore());
		}
		Map<String, Integer> a2 = new HashMap<String, Integer>();
		for (Result r : census2.getData()) {
			a2.put(r.getScopId(), r.getAlignment().getAlignLength());
		}

		// changed
		{
			double changed = 0;
			int nChanged = 0;
			for (Result r : census2.getData()) {
				if (m1.containsKey(r.getScopId())) {
					if (m1.get(r.getScopId()) != r.getAlignment().getTmScore()) {
						int nTriviallyAligned = 0;
						for (Map.Entry<Integer, Integer> entry : r.getAlignmentMapping().getSimpleFunction().entrySet()) {
							if (entry.getKey() == entry.getValue()) {
								nTriviallyAligned++;
							}
						}
						changed += r.getAlignment().getTmScore() - m1.get(r.getScopId());
						nChanged++;
						if (sig.isSignificant(r) &&  m1.get(r.getScopId()) < 0.4 && Math.abs(a1.get(r.getScopId()) - a2.get(r.getScopId())) < 0.0001) {
							System.out.println(r.getScopId() + "\t" + nTriviallyAligned + "\t" + CensusStatUtils.formatD(m1.get(r.getScopId())) + "\t" + CensusStatUtils.formatD(r.getAlignment().getTmScore()) + "\t" + a1.get(r.getScopId()) + "\t" + a2.get(r.getScopId()));
						}
					}
				}
			}
			System.out.println("Changed: " + CensusStatUtils.formatD(changed / nChanged) + "\t" + nChanged);
		}

		// added
		{
			double added = 0;
			int nAdded = 0;
			for (Result r : census2.getData()) {
				if (!m1.containsKey(r.getScopId())) {
					added += r.getAlignment().getTmScore();
					nAdded++;
				}
			}
			System.out.println("Added: " + CensusStatUtils.formatD(added / nAdded) + "\t" + nAdded);
		}

		// removed
		{
			double removed = 0;
			int nRemoved = 0;
			for (Result r : census1.getData()) {
				if (!m2.containsKey(r.getScopId())) {
					removed += r.getAlignment().getTmScore();
					nRemoved++;
				}
			}
			System.out.println("Removed: " + CensusStatUtils.formatD(removed / nRemoved) + "\t" + nRemoved);
		}
	}
	
}
