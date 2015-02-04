package org.biojava.nbio.structure.align.symm.census3.stats.order;

import java.util.Map;

public interface ConsensusDecider {

	int decide(Map<Integer,Integer> orders);
}
