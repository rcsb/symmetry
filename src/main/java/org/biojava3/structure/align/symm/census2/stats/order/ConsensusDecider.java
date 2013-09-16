package org.biojava3.structure.align.symm.census2.stats.order;

import java.util.Map;

public interface ConsensusDecider {

	int decide(Map<Integer,Integer> orders);
}
