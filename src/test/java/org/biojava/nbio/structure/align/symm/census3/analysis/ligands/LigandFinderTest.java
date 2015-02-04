package org.biojava.nbio.structure.align.symm.census3.analysis.ligands;


import junit.framework.TestCase;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.census3.analysis.ligands.LigandFinder;
import org.junit.Before;
import org.junit.Test;


/**
 * A test for {@link LigandFinder}.
 * @author dmyerstu
 */
public class LigandFinderTest extends TestCase {

	private AtomCache cache = new AtomCache();
	private static ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);

	ChemCompProvider orig;
	
	@Before
	public void setUp() throws Exception{
		super.setUp();
		
		cache.setFetchFileEvenIfObsolete(true);
		cache.setAutoFetch(true);
		
		FileParsingParameters params = cache.getFileParsingParams();
		
		params.setLoadChemCompInfo(true);
		
		params.setCreateAtomBonds(true);
		
		orig = ChemCompGroupFactory.getChemCompProvider();
		
		ChemCompProvider provider = new DownloadChemCompProvider();
		
		ChemCompGroupFactory.setChemCompProvider(provider);
		
		ScopFactory.setScopDatabase(scop);
	}
	
	@Test
	public void test() {
		
	}
	
	@Override
	protected void tearDown() throws Exception {
		super.tearDown();
		System.out.println("tear down");
		ChemCompGroupFactory.setChemCompProvider(orig);
	}
	
	@Test
	public void testNothing(){
		assertTrue(true);
	}
}
