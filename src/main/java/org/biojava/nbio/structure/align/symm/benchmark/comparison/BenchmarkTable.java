package org.biojava.nbio.structure.align.symm.benchmark.comparison;

import java.io.File;
import java.io.IOException;

import org.biojava.nbio.structure.align.symm.benchmark.KnownInfo;
import org.biojava.nbio.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusSignificance;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;

public class BenchmarkTable {

	public static void main(String[] args) throws IOException {
		if (args.length != 3) {
			System.err.println("Usage: " + BenchmarkTable.class.getSimpleName() + " ce-symm-benchmark.xml symd-published-benchmark.xml symd-unpublished-benchmark.xml");
			return;
		}
		printTable(Sample.fromXML(new File(args[0])), Sample.fromXML(new File(args[1])), Sample.fromXML(new File(args[1])));
	}

	private static void printTable(Sample cesymm, Sample symdPublished, Sample symdUnpublished) {
		CensusSignificance forCeSymmTm = CensusSignificanceFactory.forCeSymmTm();
		CensusSignificance forCeSymmOrd = CensusSignificanceFactory.forCeSymmOrd();
		CensusSignificance forSymd8 = CensusSignificanceFactory.forPublishedSymD8();
		CensusSignificance forSymd10 = CensusSignificanceFactory.forPublishedSymD10();
		CensusSignificance forSymdUnpublished = CensusSignificanceFactory.forUnpublishedSymD();
		System.out.println("----------------------");
		for (int i = 0; i < cesymm.getData().size(); i++) {
			KnownInfo info = cesymm.getData().get(i).getKnownInfo();
			CensusResult cesymmCase = cesymm.getData().get(i).getResult();
			CensusResult symdPublishedCase = symdPublished.getData().get(i).getResult();
			CensusResult symdUnpublishedCase = symdUnpublished.getData().get(i).getResult();
			System.out.print(cesymmCase.getId());
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
