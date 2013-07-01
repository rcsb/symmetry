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
 * Created on 2013-03-08
 *
 */
package org.biojava3.structure.align.symm.census2.benchmark;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.junit.Test;


/**
 * A test for {@link SampleBuilder}.
 * @author dmyerstu
 */
public class SampleBuilderTest {

	private static String RESOURCE_PATH = "src/test/resources/";
	
	@Test
	public void testBuildSample() throws IOException {
		File knownInfosFile = new File(RESOURCE_PATH + "census2/benchmark/benchmark1_known_orders");
		File inputFile = new File(RESOURCE_PATH + "census2/benchmark/benchmark1_stub.xml");
		File outputFile = new File(RESOURCE_PATH + "census2/benchmark/benchmark1_actual.xml");
		File expectedOutputFile = new File(RESOURCE_PATH + "census2/benchmark/benchmark1_expected.xml");
		Map<String,KnownInfo> knownInfos = SampleBuilder.getOrders(knownInfosFile);
		SampleBuilder.buildSample(inputFile, outputFile, knownInfos);
		Sample actual = Sample.fromXML(outputFile);
		Sample expected = Sample.fromXML(expectedOutputFile);
		assertEquals(actual.size(), 93); // note that the size is here manually
		for (int i = 0; i < actual.size(); i++) {
			if (!expected.getData().get(i).equals((actual).getData().get(i))) {
				System.err.println(actual.getData().get(i));
				fail();
			}
		}
		outputFile.delete();
	}

}
