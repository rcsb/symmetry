package org.biojava.nbio.structure.align.symm.census2;

import static org.junit.Assert.*;
import static org.junit.Assume.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.census2.Census;
import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.symm.census2.Results;
import org.biojava.nbio.structure.align.symm.census2.Significance;
import org.biojava.nbio.structure.align.symm.census2.SignificanceFactory;
import org.biojava.nbio.structure.align.symm.protodomain.ResourceList;
import org.biojava.nbio.structure.align.symm.protodomain.ResourceList.CensusResultListDifferenceListener;
import org.biojava.nbio.structure.align.symm.protodomain.ResourceList.NameProvider;
import org.custommonkey.xmlunit.DifferenceListener;
import org.junit.Before;
import org.junit.Test;

/**
 * A unit test for {@link Census}.
 * @author dmyerstu
 */
@SuppressWarnings("deprecation")
public class CensusTest {

	class TinyCensus extends Census {

		private String[] domains;
		
		public TinyCensus(String... domains) {
			super();
			this.domains = domains;
		}

		@Override
		protected Significance getSignificance() {
			return SignificanceFactory.tmScore(0.0);
		}
		
		@Override
		protected List<ScopDomain> getDomains() {
			List<ScopDomain> domains = new ArrayList<ScopDomain>();
			ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75B);
			for (String domain : this.domains) {
				domains.add(scop.getDomainByScopID(domain));
			}
			return domains;
		}

	}

	@Before
	public void setUp() throws StructureException {
		ResourceList.set(NameProvider.defaultNameProvider(), ResourceList.DEFAULT_PDB_DIR);
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75B);
		ScopFactory.setScopDatabase(scop); 
	}

	@Test
	public void testWithAlignmentMapping() throws IOException {
		File actualFile = File.createTempFile("actualresult1", "xml");
		assumeNotNull(actualFile);
		assumeTrue(actualFile.canWrite());
		Census census = new TinyCensus("d2c35e1");
		census.setCache(ResourceList.get().getCache());
		census.setOutputWriter(actualFile);
		census.setRecordAlignmentMapping(true);
		census.run();
		// unfortunately, the timestamp will be different
		DifferenceListener listener = new CensusResultListDifferenceListener();
		File expectedFile = ResourceList.get().openFile("census2/expected1_with_map.xml");
		boolean similar = ResourceList.compareXml(expectedFile, actualFile, listener);
		assertTrue(similar);
	}
	
	/**
	 * Test on live data.
	 * @throws IOException
	 */
	@Test
	public void testBasic() throws IOException {
		File actualFile = File.createTempFile("actualresult1", "xml");
		assumeNotNull(actualFile);
		assumeTrue(actualFile.canWrite());
		Census census = new TinyCensus("d2c35e1");
		census.setCache(ResourceList.get().getCache());
		census.setOutputWriter(actualFile);
		census.run();
		// unfortunately, the timestamp will be different
		DifferenceListener listener = new CensusResultListDifferenceListener();
		File expectedFile = ResourceList.get().openFile("census2/expected1.xml");
		boolean similar = ResourceList.compareXml(expectedFile, actualFile, listener);
		assertTrue(similar);
	}

	@Test
	public void testWithPartialResult() throws IOException {
		/*
		 * d2csba5
		 * d1m2oa2
		 * d1nona_
		 * d1nwub2
		 * d2p3pb_
		 */
		/*
		 * d1iwla_
		 * d1syea_
		 * d1m2vb1
		 * d1tyeb1
		 */
		String[] firstDomains = new String[] {"d2csba5", "d1m2oa2", "d1nona_", "d1nwub2", "d2p3pb_"};
		String[] secondDomains = new String[] {"d1iwla_", "d1syea_", "d1m2vb1", "d1tyeb1"};
		File actualFile = File.createTempFile("actualresult1", "xml");
		assumeNotNull(actualFile);
		assumeTrue(actualFile.canWrite());
		Census firstHalf = new TinyCensus(firstDomains);
		Census secondHalf = new TinyCensus(secondDomains);
		firstHalf.setCache(ResourceList.get().getCache());
		firstHalf.setOutputWriter(actualFile);
		secondHalf.setCache(ResourceList.get().getCache());
		secondHalf.setOutputWriter(actualFile);
		
		// run the first, then the second
		// since both have the same output file, that file should contain both
		firstHalf.run();
		secondHalf.run();
		Results results = Results.fromXML(actualFile);
		assertEquals("Wrong number of results", 5+4, results.size());
		HashSet<String> set = new HashSet<String>();
		for (String s : firstDomains) set.add(s);
		for (String s : secondDomains) set.add(s);
		for (Result result : results.getData()) {
			// this is sufficient because we know they're the same size
			assertTrue("Wrong result", set.contains(result.getScopId()));
		}
	}
	
	@Test
	public void testHard() throws IOException {
		File actualFile = File.createTempFile("actualresult2", ".xml");
		assumeNotNull(actualFile);
		assumeTrue(actualFile.canWrite());
		Census census = new TinyCensus("d1kcwa6");
		census.setCache(ResourceList.get().getCache());
		census.setOutputWriter(actualFile);
		census.run();
		// unfortunately, the timestamp will be different
		DifferenceListener listener = new CensusResultListDifferenceListener();
		File expectedFile = ResourceList.get().openFile("census2/expected2.xml");
		boolean similar = ResourceList.compareXml(expectedFile, actualFile, listener);
		assertTrue(similar);
	}
}
