/**
 * 
 */
package org.biojava3.structure.align.symm.order;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.CeSymm;
import org.junit.Before;
import org.junit.Test;

/**
 * @author spencer
 *
 */
public class RotationOrderDetectorTest {
	private RotationOrderDetector detector;
	private CeSymm ce;
	
	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {


		ce = new CeSymm();
		detector = new RotationOrderDetector(8);
		
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75A);
	}

	@Test
	public void testFitHarmonics() throws IOException, StructureException {
		String name;

		// Perform alignment to determine axis
		Atom[] ca1, ca2;
		AFPChain alignment;
		RotationAxis axis;
		double[] coefs,expectedHarmonics;

		name = "1MER.A";
		ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
		ca2 = StructureTools.cloneCAArray(ca1);
		alignment = ce.align(ca1, ca2);
		axis = new RotationAxis(alignment);

		coefs = detector.fitHarmonics(ca1, axis);
		expectedHarmonics = new double[] {
				1.218482, 2.110836, 0.7203669, 0.8226358,
				0.6092911, 0.6339138, 0.4439472, 0.4737434,
		};

		assertArrayEquals(name,expectedHarmonics,coefs,1e-4);


		ce = new CeSymm();// work around bug

		name = "d1ijqa1";
		ca1 = StructureTools.getAtomCAArray(StructureTools.getStructure(name));
		ca2 = StructureTools.cloneCAArray(ca1);
		alignment = ce.align(ca1, ca2);
		axis = new RotationAxis(alignment);

		coefs = detector.fitHarmonics(ca1, axis);
		expectedHarmonics = new double[] {
				0.5176411, 0.5359353, 0.4928912, 0.5044149,
				0.4031307, 1.915722, 0.4049375, 0.4456366
		};

		assertArrayEquals(name,expectedHarmonics,coefs,1e-4);

	}

}
