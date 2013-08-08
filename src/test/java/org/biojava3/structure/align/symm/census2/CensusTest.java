package org.biojava3.structure.align.symm.census2;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.protodomain.ResourceList;
import org.biojava3.structure.align.symm.protodomain.ResourceList.ElementTextIgnoringDifferenceListener;
import org.biojava3.structure.align.symm.protodomain.ResourceList.NameProvider;
import org.custommonkey.xmlunit.DifferenceListener;
import org.junit.Before;
import org.junit.Test;

/**
 * A unit test for {@link Census}.
 * @author dmyerstu
 */
public class CensusTest {

	class TinyCensus extends Census {

		private String[] domains;
		
		public TinyCensus(String... domains) {
			super(1); // 1 core only to ensure the results return in the expected order
			this.domains = domains;
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

	/**
	 * Test on live data.
	 * @throws IOException
	 */
	@Test
	public void testBasic() throws IOException {
		File actualFile = File.createTempFile("actualresult1", "xml");
		Census census = new TinyCensus("d2c35e1");
		census.setCache(ResourceList.get().getCache());
		census.setOutputWriter(actualFile);
		census.run();
		// unfortunately, the timestamp will be different
		DifferenceListener listener = new ElementTextIgnoringDifferenceListener("timestamp");
		File expectedFile = ResourceList.get().openFile("census2/expected1.xml");
		boolean similar = ResourceList.compareXml(expectedFile, actualFile, listener);
		assertTrue(similar);
	}

	@Test
	public void testWithPartialResult() {
		// TODO
	}
	
	@Test
	public void testHard() {
		// TODO
	}

}
