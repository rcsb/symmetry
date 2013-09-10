package org.biojava3.structure.align.symm.census2.analysis;

import static org.junit.Assert.*;

import java.util.Map;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.align.symm.census2.CensusJob;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.junit.Test;


/**
 * A test for {@link LigandFinder}.
 * @author dmyerstu
 */
public class LigandFinderTest {

	private static int RADIUS = LigandFinder.DEFAULT_RADIUS;
	
	private AtomCache cache = new AtomCache();

	@Test
	public void testInCenter1() {
		// 4-helix bundles with each heme situated in center
		String[] scopIds = new String[] {"d1hmda_", "d1hmdb_", "d1hmdc_", "d1hmdd_"};
		for (String scopId : scopIds) {
			String ligand = find(scopId, RADIUS);
			assertNotNull(ligand);
			assertEquals("*Fe2O", ligand);
		}
	}

	@Test
	public void testNotInCenter1() {
		String scopId = "d2ca5b1";
		String ligand = find(scopId, RADIUS);
		assertNull(ligand);
	}

	@Test
	public void testNotInCenter2() {
		// a 4-helix bundle with a Zinc ligand outside
		String scopId = "d1a0ba_";
		String ligand = find(scopId, RADIUS);
		assertNull(ligand);
	}

	@Test
	public void testHybrid() {
		// An organic sulfate is in the center, but a hydrocarbon is just outside
		String scopId = "d3ejba1";
		String ligand = find(scopId, RADIUS);
		assertNotNull(ligand);
		assertEquals("C8O5S", ligand); // also implies this is the only one
	}

	@Test
	public void testNoLigand() {
		// a TIM barrel with no ligand
		String scopId = "d1i45a";
		String ligand = find(scopId, RADIUS);
		assertNull(ligand);
	}
	

	@Test
	public void testAsymmetric() {
		String scopId = "d1a19a_";
		String ligand = find(scopId, RADIUS);
		assertNull(ligand);
	}
	
	private String find(String scopId, int radius) {
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
		CensusJob job = CensusJob.forScopId(Census.AlgorithmGiver.getDefault(), SignificanceFactory.rotationallySymmetricSmart(), scopId, 0, cache, scop);
		Result result = job.call();
		System.out.println(result);
		if (!SignificanceFactory.rotationallySymmetricSmart().isSignificant(result)) return null;
		Results results = new Results();
		results.add(result);
		LigandFinder finder = new LigandFinder(radius);
		finder.find(results);
		Map<String,String> formulas = finder.getFormulas();
		return formulas.get(scopId);
	}
}
