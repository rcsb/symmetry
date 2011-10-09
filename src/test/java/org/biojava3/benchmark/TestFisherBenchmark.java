package org.biojava3.benchmark;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import junit.framework.TestCase;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.client.PdbPair;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AtomCache;

public class TestFisherBenchmark extends TestCase {


	public void testBenchmark(){

		Benchmark benchmark = new FisherBenchmark();

		List<PdbPair> pairs = benchmark.getPairs();

		assertEquals(68,pairs.size());
	}

	
}
