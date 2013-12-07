package org.biojava3.structure.align.symm.benchmark.comparison;

import java.io.File;
import java.io.IOException;

import org.biojava3.structure.align.symm.benchmark.KnownInfo;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

public class BenchmarkTable {

	public static void main(String[] args) throws IOException {
		if (args.length != 3) {
			System.err.println("Usage: " + BenchmarkTable.class.getSimpleName() + " ce-symm-benchmark.xml symd-published-benchmark.xml symd-unpublished-benchmark.xml");
			return;
		}
		printTable(Sample.fromXML(new File(args[0])), Sample.fromXML(new File(args[1])), Sample.fromXML(new File(args[1])));
	}

	private static void printTable(Sample cesymm, Sample symdPublished, Sample symdUnpublished) {
		Significance forCeSymmTm = SignificanceFactory.forCeSymmTm();
		Significance forCeSymmOrd = SignificanceFactory.forCeSymmOrd();
		Significance forSymd8 = SignificanceFactory.forPublishedSymD8();
		Significance forSymd10 = SignificanceFactory.forPublishedSymD10();
		Significance forSymdUnpublished = SignificanceFactory.forUnpublishedSymD();
		System.out.println("----------------------");
		for (int i = 0; i < cesymm.getData().size(); i++) {
			KnownInfo info = cesymm.getData().get(i).getKnownInfo();
			Result cesymmCase = cesymm.getData().get(i).getResult();
			Result symdPublishedCase = symdPublished.getData().get(i).getResult();
			Result symdUnpublishedCase = symdUnpublished.getData().get(i).getResult();
			System.out.print(cesymmCase.getScopId());
			System.out.print("\t");
			System.out.print(info.getGroup());
			System.out.print("\t");
			System.out.print(forCeSymmTm.isSignificant(cesymmCase)? "T" : "F");
			System.out.print("\t");
			System.out.print(forCeSymmOrd.isSignificant(cesymmCase)? "T" : "F");
			System.out.print("\t");
			System.out.print(forSymd8.isSignificant(symdPublishedCase)? "T" : "F");
			System.out.print("\t");
			System.out.print(forSymd10.isSignificant(symdPublishedCase)? "T" : "F");
			System.out.print("\t");
			System.out.print(forSymdUnpublished.isSignificant(symdUnpublishedCase)? "T" : "F");
			System.out.println();
		}
		System.out.println("----------------------");
	}

}
