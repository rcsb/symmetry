package org.biojava3.structure.align.symm.benchmark.set;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;

import org.biojava3.structure.align.symm.benchmark.KnownInfo;
import org.biojava3.structure.align.symm.benchmark.SampleBuilder;
import org.biojava3.structure.align.symm.census3.stats.CensusStatUtils;

/**
 * Finds basic statistics about the benchmark set itself (known info).
 * @author dmyerstu
 * @since 0.2.2
 */
public class BenchmarkSetStats {

	private int nCyclic;
	private int nDihedral;
	private int nEven;
	private int nHelical;
	private int nNonIntegralHelical;
	private int nOdd;
	private int nSuperhelical;
	private int nSymm;
	private int nTotal;
	private int nTranslational;

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length > 1) {
			System.err.println("Usage: " + BenchmarkSetStats.class.getSimpleName() + " [file]");
			return;
		}
		File file = new File("src/main/resources/domain_symm_benchmark.tsv");
		if (args.length > 0) {
			file = new File(args[0]);
		}
		BenchmarkSetStats stats = new BenchmarkSetStats(file);
		System.out.println(stats);
	}

	public BenchmarkSetStats(File file) throws IOException {
		this(SampleBuilder.getOrders(file));
	}

	public BenchmarkSetStats(LinkedHashMap<String, KnownInfo> orders) {
		for (KnownInfo info : orders.values()) {
			if (info.hasDihedralSymmetry()) nDihedral++;
			if (info.hasTranslationalSymmetry()) nTranslational++;
			if (info.hasCyclicSymmetry()) nCyclic++;
			if (info.hasTrueHelicalSymmetry()) nHelical++;
			if (info.hasOddOrderSymmetry()) nOdd++;
			if (info.hasEvenOrderSymmetry()) nEven++;
			if (!info.isAsymmetric()) nSymm++;
			if (info.hasSuperhelicalSymmetry()) nSuperhelical++;
			if (info.hasNonIntegralOrderSymmetry()) nNonIntegralHelical++;
		}
		nTotal = orders.size();
	}

	public int getnCyclic() {
		return nCyclic;
	}

	public int getnDihedral() {
		return nDihedral;
	}

	public int getnEven() {
		return nEven;
	}

	public int getnHelical() {
		return nHelical;
	}

	public int getnNonIntegralHelical() {
		return nNonIntegralHelical;
	}

	public int getnOdd() {
		return nOdd;
	}

	public int getnSuperhelical() {
		return nSuperhelical;
	}

	public int getnSymm() {
		return nSymm;
	}

	public int getnTotal() {
		return nTotal;
	}

	public int getnTranslational() {
		return nTranslational;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("symmetric\t" + nSymm + CensusStatUtils.NEWLINE);
		sb.append("cyclic\t" + nCyclic + CensusStatUtils.NEWLINE);
		sb.append("dihedral\t" + nDihedral + CensusStatUtils.NEWLINE);
		sb.append("translational\t" + nTranslational + CensusStatUtils.NEWLINE);
		sb.append("true helical\t" + (nHelical - nSuperhelical - nNonIntegralHelical) + CensusStatUtils.NEWLINE);
		sb.append("superhelical\t" + nSuperhelical + CensusStatUtils.NEWLINE);
		sb.append("non-integral helical\t" + nNonIntegralHelical + CensusStatUtils.NEWLINE);
		sb.append("even-order\t" + nEven + CensusStatUtils.NEWLINE);
		sb.append("odd-order\t" + nOdd + CensusStatUtils.NEWLINE);
		sb.append("total\t" + nTotal + CensusStatUtils.NEWLINE);
		return sb.toString();
	}

}
