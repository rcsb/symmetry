/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2013-05-25
 *
 */
package org.biojava3.structure.align.symm.benchmark.comparison;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava3.structure.align.symm.benchmark.Case;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;

/**
 * Finds (significant) differences between different {@link Criterion Criteria} on either a {@link Sample} object or a {@link Results} object.
 * 
 * @author dmyerstu
 */
public class CriteriaDifferences {

	public static class Difference implements Comparable<Difference> {
		private final Double first;
		private String scopId;
		private final Double second;

		public Difference(String scopId, Double first, Double second) {
			this.scopId = scopId;
			this.first = first;
			this.second = second;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Difference other = (Difference) obj;
			if (scopId == null) {
				if (other.scopId != null)
					return false;
			} else if (!scopId.equals(other.scopId))
				return false;
			return true;
		}

		public Double getFirst() {
			return first;
		}

		public String getScopId() {
			return scopId;
		}

		public Double getSecond() {
			return second;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (scopId == null ? 0 : scopId.hashCode());
			return result;
		}

		public String print(Criterion a, Criterion b) {
			return "[" + scopId + ": " + a + "=" + first + ", " + b + "=" + second + "]";
		}

		public void setScopId(String scopId) {
			this.scopId = scopId;
		}

		@Override
		public String toString() {
			return "Difference [scopId=" + scopId + ", first=" + first + ", second=" + second + "]";
		}

		@Override
		public int compareTo(Difference o) {
			return (new Double(Math.abs(o.getFirst() - o.getSecond()))).compareTo(new Double(Math.abs(first - second)));
		}
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + CriteriaDifferences.class.getSimpleName() + " input-sample-file");
			return;
		}
		Criterion a = Criterion.tmScore();
		Criterion b = Criterion.tmpr();
		CriteriaDifferences comp = new CriteriaDifferences(a, b);
		// Results results = Results.fromXML(new File(args[0]));
		Sample sample = Sample.fromXML(new File(args[0]));
		SortedSet<Difference> differences = comp.findDifferences(sample, 0.001);
		System.out.println(differences.size() + " found:");
		for (Difference diff : differences) {
			System.out.println(diff.print(a, b));
		}
	}

	private final Criterion a;

	private final Criterion b;

	public CriteriaDifferences(Criterion a, Criterion b) {
		this.a = a;
		this.b = b;
	}

	public List<Difference> findDifferences(Results results, double precision) {
		List<Difference> diffs = new ArrayList<Difference>();
		for (Result result : results.getData()) {
			Double x = null, y = null;
			try {
				x = a.get(result);
			} catch (NoncomputableCriterionException e) {
			}
			try {
				y = b.get(result);
			} catch (NoncomputableCriterionException e) {
			}
			if (Math.abs(x) - Math.abs(y) >= precision) {
				diffs.add(new Difference(result.getScopId(), x, y));
			}
		}
		return diffs;
	}

	public SortedSet<Difference> findDifferences(Sample sample, double precision) {
		SortedSet<Difference> diffs = new TreeSet<Difference>();
		for (Case cas : sample.getData()) {
			Double x = null, y = null;
			try {
				x = a.get(cas.getResult());
			} catch (NoncomputableCriterionException e) {
				e.printStackTrace();
			}
			try {
				y = b.get(cas.getResult());
			} catch (NoncomputableCriterionException e) {
				e.printStackTrace();
			}
			if (Math.abs(x) - Math.abs(y) >= precision) {
				diffs.add(new Difference(cas.getScopId(), x, y));
			}
		}
		return diffs;
	}

}
