package org.biojava3.structure.align.symm.census2.benchmark;


import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.protodomain.ResourceList;
import org.biojava3.structure.align.symm.protodomain.ResourceList.NameProvider;
import org.junit.Before;
import org.junit.Test;

public class SampleTest {

	@Before
	public void setUp() throws Exception {
		ResourceList.set(NameProvider.defaultNameProvider(), ResourceList.DEFAULT_PDB_DIR);
		ScopDatabase scop = ScopFactory.getSCOP();
		if (!scop.getClass().getName().equals(BerkeleyScopInstallation.class.getName())) { // for efficiency
			ScopFactory.setScopDatabase(new BerkeleyScopInstallation()); // ScopDatabase is too hard to mock well
		}
	}

	@Test
	public void testRead() throws IOException {
		Sample sample = Sample.fromXML(ResourceList.get().openFile("census2/benchmark/test_read.xml"));
		List<Case> cases = sample.getData();
		Case c = cases.get(0);
		assertEquals(1, cases.size());
		assertEquals("d1ezda1", c.getScopId());
		assertEquals("D4", c.getKnownGroup());
		assertEquals(4, c.getKnownOrder());
		assertEquals("1ezd.A_40-62,A_68-144", c.getProtodomain());
		assertEquals(4.74, c.getAlignment().getzScore(), 0.000001);
		assertEquals(0.64949524, c.getAxis().getScrew(), 0.000001);
	}

}
