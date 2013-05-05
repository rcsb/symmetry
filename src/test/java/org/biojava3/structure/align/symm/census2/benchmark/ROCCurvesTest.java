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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava3.structure.align.symm.census2.Result;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

/**
 * A test for {@link ROCCurves}.
 * @author dmyerstu
 */
public class ROCCurvesTest {

	private static String RESOURCE_PATH = "src/test/resources/";
	
	@Test
	public void testGetRocPoints() throws IOException {
		
		File inputFile = new File(RESOURCE_PATH + "census2/benchmark/benchmark2.xml");
		List<Criterion> criteria = new ArrayList<Criterion>();
		criteria.add(Criterion.alignLength());
		criteria.add(Criterion.screw());
		criteria.add(new Criterion() {
			@Override
			public double get(Result result) throws NoncomputableCriterionException {
				if (result.getAlignment() == null || result.getAlignment().getInitialShift() == null) throw new NoncomputableCriterionException();
				return result.getAlignment().getInitialShift();
			}
			@Override
			public String getName() {
				return "initial shift";
			}
		});
		
		ROCCurves rocs = new ROCCurves(inputFile, criteria);
		XYSeriesCollection dataset = rocs.getRocPoints();
		
		// we shouldn't have more criteria
		assertEquals(criteria.size(), dataset.getSeriesCount());
		
		float numPos, numNeg;
		
		/*
		 * Aligned length:
		 * T=90
		 * T=80
		 * F=70
		 * T=60
		 * F=40
		 * F=30
		 * T=20
		 */
		numPos = 3.0f;
		 numNeg = 4.0f;
		XYSeries alignLength = dataset.getSeries(0);
		assertEquals(7, alignLength.getItemCount());
		
		assertEquals(1, alignLength.getY(0).floatValue()*numPos, 0.001f);
		assertEquals(0, alignLength.getX(0).floatValue()*numNeg, 0.001f);
		
		assertEquals(2, alignLength.getY(1).floatValue()*numPos, 0.001f);
		assertEquals(0, alignLength.getX(1).floatValue()*numNeg, 0.001f);

		assertEquals(2, alignLength.getY(2).floatValue()*numPos, 0.001f);
		assertEquals(1, alignLength.getX(2).floatValue()*numNeg, 0.001f);
		
		assertEquals(3, alignLength.getY(3).floatValue()*numPos, 0.001f);
		assertEquals(1, alignLength.getX(3).floatValue()*numNeg, 0.001f);
		
		assertEquals(3, alignLength.getY(4).floatValue()*numPos, 0.001f);
		assertEquals(2, alignLength.getX(4).floatValue()*numNeg, 0.001f);

		assertEquals(3, alignLength.getY(5).floatValue()*numPos, 0.001f);
		assertEquals(3, alignLength.getX(5).floatValue()*numNeg, 0.001f);

		assertEquals(3, alignLength.getY(6).floatValue()*numPos, 0.001f);
		assertEquals(4, alignLength.getX(6).floatValue()*numNeg, 0.001f);
		

		/*
		 * Screw:
		 * T=11
		 * T=5
		 * F=1
		 * F=-10
		 */
		XYSeries screw = dataset.getSeries(1);
		assertEquals(4, screw.getItemCount());
		numPos = 2.0f;
		numNeg = 2.0f;

		assertEquals(1, screw.getY(0).floatValue()*numPos, 0.001f);
		assertEquals(0, screw.getX(0).floatValue()*numNeg, 0.001f);

		assertEquals(2, screw.getY(1).floatValue()*numPos, 0.001f);
		assertEquals(0, screw.getX(1).floatValue()*numNeg, 0.001f);

		assertEquals(2, screw.getY(2).floatValue()*numPos, 0.001f);
		assertEquals(1, screw.getX(2).floatValue()*numNeg, 0.001f);
		
		assertEquals(2, screw.getY(3).floatValue()*numPos, 0.001f);
		assertEquals(2, screw.getX(3).floatValue()*numNeg, 0.001f);

		/*
		 * Screw:
		 * F=20
		 * T=-50
		 * T=-80
		 */
		XYSeries shift = dataset.getSeries(2);
		assertEquals(3, shift.getItemCount());
		numPos = 2.0f;
		numNeg = 1.0f;

		assertEquals(0, shift.getY(0).floatValue()*numPos, 0.001f);
		assertEquals(1, shift.getX(0).floatValue()*numNeg, 0.001f);

		assertEquals(1, shift.getY(1).floatValue()*numPos, 0.001f);
		assertEquals(1, shift.getX(1).floatValue()*numNeg, 0.001f);

		assertEquals(2, shift.getY(2).floatValue()*numPos, 0.001f);
		assertEquals(1, shift.getX(2).floatValue()*numNeg, 0.001f);
		
	}
	
}
