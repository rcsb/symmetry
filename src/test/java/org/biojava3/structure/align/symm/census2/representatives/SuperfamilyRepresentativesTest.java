package org.biojava3.structure.align.symm.census2.representatives;


import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.protodomain.ResourceList;
import org.biojava3.structure.align.symm.protodomain.ResourceList.NameProvider;
import org.junit.Before;
import org.junit.Test;


public class SuperfamilyRepresentativesTest {

	@Before
	public void setUp() throws StructureException {
		ResourceList.set(NameProvider.defaultNameProvider(), ResourceList.DEFAULT_PDB_DIR);
		ScopFactory.setScopDatabase(ScopFactory.getSCOP(ScopFactory.VERSION_1_75A)); // the test will break for a different SCOP version
	}

	@Test
	public void testWithFamily() {
		// b.1.1.3 and b.1.1.4
		SuperfamilyRepresentatives reps = new SuperfamilyRepresentatives(3, new int[] {49159, 49142});
		List<ScopDomain> domains = reps.getDomains();
		assertEquals(6, domains.size());
		assertEquals("b.1.1.4", domains.get(0).getClassificationId());
		assertEquals("b.1.1.4", domains.get(1).getClassificationId());
		assertEquals("b.1.1.4", domains.get(2).getClassificationId());
		assertEquals("b.1.1.3", domains.get(3).getClassificationId());
		assertEquals("b.1.1.3", domains.get(4).getClassificationId());
		assertEquals("b.1.1.3", domains.get(5).getClassificationId());
	}

	@Test
	public void testWithDomain() {
		SuperfamilyRepresentatives reps = new SuperfamilyRepresentatives(100, new int[] {69160, 49196});
		List<ScopDomain> domains = reps.getDomains();
		assertEquals(2, domains.size()); // since we have only 2 domains, we can't find 100
	}

	@Test
	public void testWithSpecies() {
		SuperfamilyRepresentatives reps = new SuperfamilyRepresentatives(null, new int[] {49183});
		List<ScopDomain> domains = reps.getDomains();
		assertEquals(1, domains.size());
		// the first Px is 21752
		assertEquals(21752, domains.get(0).getSunid().intValue());
	}

	@Test
	public void testWithPx() {
		SuperfamilyRepresentatives reps = new SuperfamilyRepresentatives(null, new int[] {21752});
		List<ScopDomain> domains = reps.getDomains();
		assertEquals(1, domains.size());
		assertEquals(21752, domains.get(0).getSunid().intValue());
	}

	/**
	 * <strong>Warning: this test takes a very long time</strong>.
	 */
	@Test
	public void verify1PerFamily() {
		SuperfamilyRepresentatives reps = new SuperfamilyRepresentatives();
		List<ScopDomain> domains = reps.getDomains();
		assertEquals(10270, domains.size());
		HashMap<String,ScopDescription> sfs = reps.getSuperfamilies();
		HashSet<Integer> s = new HashSet<Integer>();
		for (ScopDescription d : sfs.values()) s.add(d.getSunID());
		assertEquals(1961-129, s.size());
	}

	/**
	 * <strong>Warning: this test takes a very long time</strong>.
	 */
	@Test
	public void verify1PerSf() {
		SuperfamilyRepresentatives reps = new SuperfamilyRepresentatives(1);
		List<ScopDomain> domains = reps.getDomains();
		assertEquals(1832, domains.size());
		Set<Integer> includedSfs = new HashSet<Integer>();
		// we should never have more than 1 of the same
		for (ScopDomain domain : domains) {
			int sf = domain.getSuperfamilyId();
			assertFalse(includedSfs.contains(sf));
			includedSfs.add(sf);
		}
		HashMap<String,ScopDescription> sfs = reps.getSuperfamilies();
		assertEquals(1832, sfs.size());
		for (ScopDescription d : sfs.values()) {
			assertTrue(includedSfs.contains(d.getSunID()));
		}
	}

	/**
	 * <strong>Warning: this test takes a very long time</strong>.
	 */
	@Test
	public void verify2PerSf() {
		SuperfamilyRepresentatives reps = new SuperfamilyRepresentatives(2);
		List<ScopDomain> domains = reps.getDomains();
		assertEquals(2917, domains.size());
		Map<Integer,Integer> count = new HashMap<Integer,Integer>();
		Set<Integer> includedFamilies = new HashSet<Integer>();
		for (ScopDomain domain : domains) {
			int sf = domain.getFamilyId();
			if (!count.containsKey(sf)) count.put(sf, 0);
			includedFamilies.add(sf);
			count.put(sf, count.get(sf) + 1); // increment count
		}
		// test count is always 2
		for (int n : count.values()) {
			assertTrue("Wrong number of domains", n == 1 || n == 2);
		}
		
	}
	
}
