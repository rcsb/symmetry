/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2013-03-01
 *
 */
package org.biojava.nbio.structure.align.symm.benchmark;


import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.scop.ScopFactory;
import org.junit.Before;
import org.junit.Test;

/**
 * A test for {@link Sample}.
 * @author dmyerstu
 */
public class SampleTest {

	private static String RESOURCE_PATH = "src/test/resources/";
	
	@Before
	public void setUp() throws Exception {
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75A);
	}

	@Test
	public void testRead() throws IOException {
		Sample sample = Sample.fromXML(new File(RESOURCE_PATH + "census2/benchmark/test_read.xml"));
		List<Case> cases = sample.getData();
		Case c = cases.get(0);
		assertEquals(1, cases.size());
		assertEquals("d1ezda1", c.getScopId());
		assertEquals("D4", c.getKnownGroup());
		assertEquals(4, c.getKnownOrder());
		assertEquals("1ezd.A_40-62,A_68-144", c.getAlignedUnit());
		assertEquals(4.74, c.getScoreList().getzScore(), 0.000001);
		assertEquals(0.64949524, c.getAxis().getParallel(), 0.000001);
	}

}
