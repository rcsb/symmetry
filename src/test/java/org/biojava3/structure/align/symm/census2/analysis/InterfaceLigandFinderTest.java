package org.biojava3.structure.align.symm.census2.analysis;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.align.symm.census2.CensusJob;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.junit.Test;

/**
 * A test for {@link InterfaceLigandFinder}.
 * @author dmyersturnbull
 * @deprecated
 */
@Deprecated
public class InterfaceLigandFinderTest {

	private AtomCache cache = new AtomCache();
	private ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75B);

	@Test
	public void nearInCenterOfAxis() {
		String name = "1hmd.A";
		String ligand = find(name, 5, false, SignificanceFactory.ultraLiberal());
		assertNotNull(ligand);
	}

	private String find(String name, int radius, boolean useOnlyAligned, Significance sig) {
		Result result = CensusJob.runJob(name, 0, Census.AlgorithmGiver.getDefault(), SignificanceFactory.rotationallySymmetricSmart(), cache, scop);
		System.out.println(result);
		Results results = new Results();
		results.add(result);
		InterfaceLigandFinder finder = new InterfaceLigandFinder();
		finder.setMaxDistance(radius);
		finder.setSignificance(sig);
		finder.find(results);
		Map<String,String> formulas = finder.getFormulas();
		if (("[" + name + "]").equals(formulas.get(name))) return null;
		return formulas.get(name);
	}
}
