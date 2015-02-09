package org.biojava.nbio.benchmark;


import java.util.List;

import junit.framework.TestCase;

import org.biojava.nbio.structure.align.client.PdbPair;
import org.biojava.nbio.benchmark.Benchmark;
import org.biojava.nbio.benchmark.FisherBenchmark;


public class TestFisherBenchmark extends TestCase {


	public void testBenchmark(){

		Benchmark benchmark = new FisherBenchmark();

		List<PdbPair> pairs = benchmark.getPairs();

		assertEquals(68,pairs.size());
	}

	
}
