package org.biojava.nbio.benchmark;

import java.util.List;

import org.biojava.nbio.structure.align.client.PdbPair;

public interface Benchmark {
	
	public String getName();
	
	public List<PdbPair> getPairs();

}
