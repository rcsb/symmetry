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
 * Created on 2013-03-22
 *
 */
package org.biojava.nbio.structure.align.symm.census2;

import org.biojava.nbio.structure.align.symm.census2.SignificanceFactory;
import org.junit.Test;

/**
 * A simple test for the reflection-based methods in {@link SignificanceFactory}.
 * @author dmyerstu
 */
public class SignificanceFactoryTest {

	@Test
	public void testFromMethod() {
		SignificanceFactory.fromMethod(null, "liberalTmScore");
	}

	@Test
	public void testFromMethodWithArgs() {
		SignificanceFactory.fromMethod(null, "tmScore", new Object[]{0.3});
	}

	@Test
	public void testFromMethodAndClass() {
		SignificanceFactory.fromMethod("org.biojava.nbio.structure.align.symm.census2.SignificanceFactory", "liberalTmScore", new Object[]{});
	}

	@Test
	public void testFromClassNoArgs() {
		SignificanceFactory.fromClass("org.biojava.nbio.structure.align.symm.census2.FakeSignificance", new Object[]{});
	}

	@Test
	public void testFromClassWithArgs() {
		SignificanceFactory.fromClass("org.biojava.nbio.structure.align.symm.census2.FakeSignificance", new Object[]{"thisisastring"});
	}

	@Test(expected=IllegalArgumentException.class)
	public void testFromClassWithArgsFail() {
		SignificanceFactory.fromClass("org.biojava.nbio.structure.align.symm.census2.FakeSignificance", new Object[]{1});
	}
	
}
