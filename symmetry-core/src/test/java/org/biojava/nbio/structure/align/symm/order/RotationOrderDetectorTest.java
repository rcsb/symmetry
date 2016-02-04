/**
 * 
 */
package org.biojava.nbio.structure.align.symm.order;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;

import static org.biojava.nbio.structure.align.symm.order.RotationOrderDetector.RotationOrderMethod.*;

import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.RefinerFailedException;
import org.junit.Before;
import org.junit.Test;

/**
 * Test all the methods for order detection.
 * 
 * @author Spencer Bliven
 * 
 */
public class RotationOrderDetectorTest {

	private CESymmParameters params;

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {

		params = new CESymmParameters();
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75A);
	}

	@Test
	public void testFitHarmonics() throws IOException, StructureException {
		String name;

		RotationOrderDetector detector = new RotationOrderDetector(8, HARMONICS);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs, expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.tryAllOrders(ca1, axis, false);
		expectedHarmonics = new double[] { 0, 1.218482, 2.110836, 0.7203669,
				0.8226358, 0.6092911, 0.6339138, 0.4439472, 0.4737434, };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

		name = "d1ijqa1";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.tryAllOrders(ca1, axis, false);
		expectedHarmonics = new double[] { 0, 0.5176411, 0.5359353, 0.4928912,
				0.5044149, 0.4031307, 1.915722, 0.4049375, 0.4456366 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

	}

	@Test
	public void testCalculateOrderByHarmonics() throws IOException,
			StructureException, RefinerFailedException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8, HARMONICS);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		int order;

		name = "1MER.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		order = detector.calculateOrder(alignment, ca1);
		assertEquals(name, 2, order);

		name = "d1ijqa1";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		order = detector.calculateOrder(alignment, ca1);
		assertEquals(name, 6, order);

		name = "1TIM.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		order = detector.calculateOrder(alignment, ca1);
		assertEquals(name, 1, order);// tough case
	}

	@Test
	public void testFitHarmonicsFloating() throws IOException,
			StructureException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8,
				HARMONICS_FLOATING);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs, expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.tryAllOrders(ca1, axis, true);
		expectedHarmonics = new double[] { 2.597383371368633,
				0.5213221707509845, 1.4121594197036587, 0.06167606937029615,
				0.17901478726304354, 0.019300560598982656, 0.06908612977051289,
				-0.06104875791630157, -0.003831898183546225 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

		name = "d1ijqa1";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.tryAllOrders(ca1, axis, true);
		expectedHarmonics = new double[] { 1.2739427717435365,
				0.2761280960318058, 0.2718581473051768, 0.22837118857649835,
				0.20840059641118155, 0.10054298829557395, 1.5787950609554844,
				0.06204472152219612, 0.07230120934912426 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

		name = "1TIM.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.tryAllOrders(ca1, axis, true);
		expectedHarmonics = new double[] { 2.835403848327179,
				0.1468762702349724, 0.1159408606767334, 0.09051681783286569,
				0.10420798340774699, -0.04381730637594528, 0.2023844545733075,
				-0.0059047948160164, 0.17631879026416564 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

	}

	@Test
	public void testCalculateOrderByHarmonicsFloating() throws IOException,
			StructureException, RefinerFailedException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8,
				HARMONICS_FLOATING);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		int order;

		name = "1MER.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		order = detector.calculateOrder(alignment, ca1);
		assertEquals(name, 2, order);

		name = "d1ijqa1";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		order = detector.calculateOrder(alignment, ca1);
		assertEquals(name, 6, order);

		name = "1TIM.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		order = detector.calculateOrder(alignment, ca1);
		assertEquals(name, 6, order);// tough case

	}

	@Test
	public void testTrySingleHarmonicsFloatingByAmp() throws IOException,
			StructureException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8,
				SINGLE_HARMONIC_AMP);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs, expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersByAmp(ca1, axis);
		expectedHarmonics = new double[] { 0.009384193201414186,
				1.1740333517311115, -0.44686757545089734, -0.2884995929083228,
				-0.2677652242044436, -0.19546335129315567,
				-0.15476867112456613, -0.09201871044036189 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

		name = "d1ijqa1";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersByAmp(ca1, axis);
		expectedHarmonics = new double[] { -0.09799681735927843,
				-0.21364824515873054, -0.11035897329374683,
				-0.21342288908870446, -0.1855205123692191, 1.4483981845191647,
				-0.14753158983373432, -0.20318262394650666 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

		name = "1TIM.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersByAmp(ca1, axis);
		expectedHarmonics = new double[] { 0.05425570731331745,
				-0.019746339856520342, -0.007014170110096962,
				-0.010253935127542774, -0.1363341693467633,
				0.13180804348846586, -0.06217377718239171, 0.13580312417478568 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

	}

	@Test
	public void testTrySingleHarmonicsFloatingBySSE() throws IOException,
			StructureException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8,
				SINGLE_HARMONIC_SSE);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs, expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersBySSE(ca1, axis);
		expectedHarmonics = new double[] { 0.4072455479103931,
				0.14067952827982602, 0.3779933525800668, 0.39480256016817444,
				0.3960406707063868, 0.4012976706684556, 0.4035501028422771,
				0.4059768120251738 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

	}

	@Test
	public void testTrySingleCuspByAmp() throws IOException, StructureException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8,
				SINGLE_CUSP_AMP);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs, expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersByAmp(ca1, axis);
		expectedHarmonics = new double[] { 0.1820283077630485,
				1.0136815511756285, -0.3905901172609695, -0.26201199397108704,
				-0.20690437223329403, -0.1616881748238757,
				-0.11567568721588678, -0.08632416298270444 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

		name = "d1ijqa1";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersByAmp(ca1, axis);
		expectedHarmonics = new double[] { -0.1476801436899504,
				-0.12885346994573413, 0.16933007895862143,
				-0.21743126249158864, -0.18435161092570687, 1.1878594622938643,
				-0.13766536861375353, -0.1779357492984687 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

		name = "1TIM.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersByAmp(ca1, axis);
		expectedHarmonics = new double[] { 0.048573997778241465,
				-0.004122183485246968, 0.016836424264462534,
				0.011565079482636204, -0.11506791209963799,
				0.10046269754081401, -0.055177206937627045, 0.10338290890177872 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

	}

	@Test
	public void testTrySingleCuspBySSE() throws IOException, StructureException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8,
				SINGLE_CUSP_SSE);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs, expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersBySSE(ca1, axis);
		expectedHarmonics = new double[] { 0.4023477, 0.1620938, 0.3740037,
				0.3923604, 0.3973309, 0.4011318, 0.4041583, 0.4054583 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);

	}

	@Test
	public void testTrySingleCuspFixedByAmp() throws IOException,
			StructureException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8,
				SINGLE_CUSP_FIXED_AMP);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs, expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersByAmp(ca1, axis);
		expectedHarmonics = new double[] { 0.1287134, 1.030371, -0.5324238,
				-0.4618442, -0.5114463, -0.4755183, -0.4395435, -0.3174656 };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-4);
	}

	@Test
	public void testTrySingleCuspFixedBySSE() throws IOException,
			StructureException {
		String name;
		RotationOrderDetector detector = new RotationOrderDetector(8,
				SINGLE_CUSP_FIXED_SSE);

		// Perform alignment to determine axis
		Atom[] ca1;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs, expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
		alignment = CeSymm.analyze(ca1, params).getSelfAlignment();
		axis = new RotationAxis(alignment);

		coefs = detector.trySingleOrdersBySSE(ca1, axis);
		expectedHarmonics = new double[] { 0.4023477, 0.1485084, 0.3794772,
				0.3946517, 0.3969921, 0.4006527, 0.4033445, 0.4056923, };

		assertArrayEquals(name, expectedHarmonics, coefs, 1e-2);
	}

}
