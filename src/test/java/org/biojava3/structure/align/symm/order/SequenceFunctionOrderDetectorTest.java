package org.biojava3.structure.align.symm.order;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.align.symm.CeSymm;
import org.junit.Test;

/**
 * Originally part of {@link CeSymmTest}.
 * @author Spencer Bliven
 */
public class SequenceFunctionOrderDetectorTest {

	@Test
	public void testGetSymmetryOrder() throws IOException, StructureException, OrderDetectionFailedException {
		// List of alignments to try, along with proper symmetry
		Map<String,Integer> orderMap = new HashMap<String,Integer>();
		orderMap.put("1itb.A",3); // b-trefoil, C3
		orderMap.put("1tim.A",2); // tim-barrel, C8
		//orderMap.put("d1p9ha_",-1); // not rotational symmetry
//		orderMap.put("3HKE.A",2); // very questionable alignment
		orderMap.put("d1jlya1",3); // a very nice trefoil
		
		AtomCache cache = new AtomCache();
		
		for(String name : orderMap.keySet()) {
			CeSymm ce = new CeSymm();
			
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			AFPChain afpChain = ce.align(ca1, ca2);
			
			int order = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
			
			assertEquals("Wrong order for "+name,orderMap.get(name).intValue(), order);
		}
	}
	
}
