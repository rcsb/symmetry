package org.biojava3.structure.align.symm.census2.analysis;


import java.util.Map;

import junit.framework.TestCase;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.ChemCompProvider;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.align.symm.census2.CensusJob;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.junit.Before;
import org.junit.Test;


/**
 * A test for {@link LigandFinder}.
 * @author dmyerstu
 */
public class LigandFinderTest extends TestCase{

	private static int RADIUS = 5;
	
	private AtomCache cache = new AtomCache();
	private ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75B);

	ChemCompProvider orig ;
	
	@Before
	public void setUp() throws Exception{
		super.setUp();
		
		cache.setFetchFileEvenIfObsolete(true);
		cache.setAutoFetch(true);
		
		FileParsingParameters params = cache.getFileParsingParams();
		
		params.setLoadChemCompInfo(true);
		
		params.setCreateAtomBonds(true);
		
		orig =  ChemCompGroupFactory.getChemCompProvider();
		
		ChemCompProvider provider = new DownloadChemCompProvider();
		
		ChemCompGroupFactory.setChemCompProvider(provider);
		
		
		
	}
	
	

	@Override
	protected void tearDown() throws Exception {
		// TODO Auto-generated method stub
		super.tearDown();
		
		System.out.println("tear down");
		ChemCompGroupFactory.setChemCompProvider(orig);
	}
	
	



	@Test
	public void testInCenter1() {
		// 4-helix bundles with each heme situated in center
		String[] scopIds = new String[] {"d1hmda_", "d1hmdb_", "d1hmdc_", "d1hmdd_"};
		for (String scopId : scopIds) {
			String ligand = find(scopId, RADIUS, false, SignificanceFactory.forCeSymmOrd());
			assertNotNull("Found no ligand for " + scopId, ligand);
			assertTrue(ligand.contains("Fe2O"));
		}
	}

	@Test
	public void testNotNearAxis1() {
		// ligand is near centroid but not near axis
		String name = "2VR1.B"; // could try 4BCT
		String centroid = find(name, 5, false, SignificanceFactory.ultraLiberal());
		assertNotNull(centroid);
		String symmetry = find(name, 5, true, SignificanceFactory.ultraLiberal());
		assertNull(symmetry);
	}
	
	@Test
	public void testNearAndFar() {
		// A transcriptional regulator (from paper)
		String scopId = "d3ddva1";
		String ligand = find(scopId, 10, false, SignificanceFactory.forCeSymmOrd());
		assertNotNull(ligand);
		assertFalse(ligand.contains(","));
		assertTrue(ligand.contains("Mg"));
		ligand = find(scopId, RADIUS, false, SignificanceFactory.forCeSymmOrd());
		assertNull(ligand);
	}

	@Test
	public void testNotInCenter1() {
		// callogen; ligands are actually pretty close to centroid
		String scopId = "d1caga_";
		String ligand = find(scopId, RADIUS, false, SignificanceFactory.forCeSymmOrd());
		assertNull(ligand);
	}

	@Test
	public void testNotInCenter2() {
		// a 4-helix bundle with a Zinc ligand outside
		String scopId = "d1a0ba_";
		String ligand = find(scopId, RADIUS, false, SignificanceFactory.forCeSymmOrd());
		assertNull(ligand);
	}

	@Test
	public void testNotInCenter3() {
		// symmetry along interface (from paper)
		String scopId = "d1squa_";
		String ligand = find(scopId, RADIUS, false, SignificanceFactory.forCeSymmOrd());
		assertNull(ligand);
	}

	@Test
	public void testHybrid() {
		// An organic sulfate is in the center, but a hydrocarbon is just outside
		String scopId = "d3ejba1";
		String ligand = find(scopId, RADIUS, false, SignificanceFactory.forCeSymmOrd());
		assertNotNull(ligand);
		assertFalse(ligand.contains(","));
		assertTrue(ligand.contains("C8O5S"));
	}

	@Test
	public void testHybridOnlyAligned() {
		String scopId = "d3ejba1";
		String ligand = find(scopId, RADIUS, true, SignificanceFactory.forCeSymmOrd());
		assertNull(ligand);
	}

	@Test
	public void testNoLigand() {
		// a TIM barrel with no ligand
		String scopId = "d1ypia_";
		String ligand = find(scopId, RADIUS, false, SignificanceFactory.forCeSymmOrd());
		assertNull(ligand);
	}

	@Test
	public void testAsymmetric() {
		String scopId = "d1hmda_";
		String ligand = find(scopId, RADIUS, false, SignificanceFactory.not(SignificanceFactory.ultraLiberal()));
		assertNull(ligand);
	}
	
	private String find(String name, int radius, boolean useOnlyAligned, Significance sig) {
		Result result = CensusJob.runJob(name, 0, Census.AlgorithmGiver.getDefault(), SignificanceFactory.rotationallySymmetricSmart(), cache, scop);
		System.out.println(result);
		Results results = new Results();
		results.add(result);
		LigandFinder finder = new LigandFinder(radius);
		finder.setUseOnlyAligned(useOnlyAligned);
		finder.setSignificance(sig);
		finder.find(results);
		Map<String,String> formulas = finder.getFormulas();
		if (("[" + name + "]").equals(formulas.get(name))) return null;
		return formulas.get(name);
	}
}
