package org.biojava3.structure.align.symm.census2.analysis.mespeus;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.census2.analysis.mespeus.SymmetryInMespeus.MespeusEntryMatcher;
import org.junit.Test;

/**
 * An integration test for {@link SymmetryInMespeus}.
 * @author dmyersturnbull
 */
public class SymmetryInMespeusTest {

	@Test
	public void testRead() throws IOException {
		SymmetryInMespeus mespeus = new SymmetryInMespeus(new File("src/test/resources/mespeus_ex.tsv"), SignificanceFactory.ultraLiberal());
		MespeusEntry first = mespeus.getEntries().get(0);
		assertEquals("1be7", first.getPdbId());
		assertEquals(4, first.getCoordinationNumber());
		assertEquals(1.937, first.getDistance(), 0.0000001f);
		CoordinationGeometry geometry = first.getShape();
		int i = 0;
		CoordinationGeometryType[] types = new CoordinationGeometryType[] {CoordinationGeometryType.TETRAHEDRAL, CoordinationGeometryType.SQUARE_PLANAR};
		float[] deltas = new float[] {3.0f, 42.4f};
		for (Map.Entry<CoordinationGeometryType, Float> entry : geometry.getDeltas().entrySet()) {
			assertEquals(types[i], entry.getKey());
			assertEquals(deltas[i], (float) entry.getValue(), 0.0000001f);
			i++;
		}
		for (MespeusEntry entry : mespeus.getEntries()) {
			System.out.println(entry);
		}
		DescriptiveStatistics stats = mespeus.correlate(Results.fromXML("src/test/resources/mespeus_census.xml"), MespeusEntryMatcher.everything());
		assertEquals(9, stats.getN());
		assertEquals(5.0 / 9.0, stats.getMean(), 0.0000001);
	}
	
}
