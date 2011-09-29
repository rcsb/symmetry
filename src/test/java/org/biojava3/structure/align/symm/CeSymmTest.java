/**
 * 
 */
package org.biojava3.structure.align.symm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;

/**
 * @author Spencer Bliven
 *
 */
public class CeSymmTest extends TestCase {

	/* (non-Javadoc)
	 * @see junit.framework.TestCase#setUp()
	 */
	protected void setUp() throws Exception {
		super.setUp();
	}

	public void testGetSymmetryOrder() throws IOException, StructureException {
		// List of alignments to try, along with proper symmetry
		Map<String,Integer> orderMap = new HashMap<String,Integer>();
		orderMap.put("1itb.A",3); // b-trefoil, C3
		orderMap.put("1tim.A",2); // tim-barrel, C8
		orderMap.put("d1p9ha_",-1); // not rotational symmetry
		orderMap.put("3HKE.A",2); // very questionable alignment
		orderMap.put("d1jlya1",3); // a very nice trefoil
		
		AtomCache cache = new AtomCache();
		//StructureAlignmentFactory.addAlgorithm(new CeSymm());
		//CeSymm ce = (CeSymm) StructureAlignmentFactory.getAlgorithm(CeSymm.algorithmName);
		
		for(String name : orderMap.keySet()) {
			CeSymm ce = new CeSymm();
			
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			AFPChain afpChain = ce.align(ca1, ca2);
			
			int order = CeSymm.getSymmetryOrder(afpChain);
			
			assertEquals("Wrong order for "+name,orderMap.get(name).intValue(), order);
		}
	}
	
	
	/**
	 * CeSymm and CeCalculator both use some internal variables, which should
	 * be reset at each call of align. Test that the order of align calls doesn't
	 * impact the results
	 * @throws StructureException 
	 * @throws IOException 
	 */
	
	/* TODO Known to fail. Fix the bug, then re-enable the test
	public void testAlignOrderIndependence() throws IOException, StructureException {
		String[] names = new String[] {
				"1itb.A", // b-trefoil, C3
				//"1tim.A", // tim-barrel, C8
				//"d1p9ha_", // not rotational symmetry
				//"3HKE.A", // very questionable alignment
				"d1jlya1", // a very nice trefoil
		};
		
		CeSymm ce;
		AtomCache cache = new AtomCache();
		AFPChain[] independent = new AFPChain[names.length];
		
		for(int i=0;i<names.length;i++) {
			// clean ce instance
			ce = new CeSymm();
			
			Atom[] ca1 = cache.getAtoms(names[i]);
			Atom[] ca2 = cache.getAtoms(names[i]);
			
			independent[i] = ce.align(ca1, ca2);
		}
		
		// Some different orders to test
		List<List<Integer>> orders = new ArrayList<List<Integer>>();
		List<Integer> inOrder = new ArrayList<Integer>();
		List<Integer> reverseOrder = new ArrayList<Integer>();
		for(int i=0;i<names.length;i++) {
			//in-order
			inOrder.add(i);
			reverseOrder.add(names.length-i-1);
		}
		orders.add(inOrder);
		orders.add(reverseOrder);
		
		// Try some random orders
		final int numRandom = 0;
		for(int i=0;i<numRandom;i++) {
			List<Integer> random = new ArrayList<Integer>(inOrder);
			Random rand = new Random(1234l); //constant seed for reproducability
			Collections.shuffle(random,rand);
			orders.add(random);
		}
		
		//Test all orders
		for(List<Integer> order : orders) {
			//Single CE per order
			ce = new CeSymm();
			
			for(int i : order) {
				String name = names[i];
				
				Atom[] ca1 = cache.getAtoms(names[i]);
				Atom[] ca2 = cache.getAtoms(names[i]);

				AFPChain afpChain = ce.align(ca1, ca2);
				
				assertEquals("Found a different result in protein "+i+" in order "+order,
						independent[i], afpChain);
			}		
		}
	}
	
	*/
}
