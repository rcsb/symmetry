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
package org.biojava3.structure;

import static org.junit.Assert.*;

import java.util.List;

import org.biojava.bio.structure.scop.ScopDomain;
import org.junit.Test;


public class MockScopDatabaseTest {

//	@Test
	public void testSimple() {
		MockScopDatabase scop = new MockScopDatabase();
		scop.addClaLine("d2c35e1	2c35	E:14-139	a.60.8.2	129717	cl=46456,cf=47768,sf=47819,fa=69044,dm=69045,sp=140642,px=129717");
		scop.addDesLine("129717	px	a.60.8.2	d2c35e1	2c35 E:14-139");
		scop.addHieLine("140642	69045	129717"); // 129711,129714,129720
		scop.addHieLine("69045	69044	140642"); // 116939,69046
		scop.addHieLine("69044	47819	69045");
		scop.addHieLine("47819	47768	69044"); // 47820,140643,140646
		scop.addHieLine("47768	46456	47819"); // 47769,47781,47789,47794,47798,47802,47807,47823,47831,69047,81585,81799,116742,140652,158544
		scop.addHieLine("46456	0	47768");
		scop.addHieLine("0	-	46456");
		List<ScopDomain> fromSunId = scop.getScopDomainsBySunid(129717);
		assertEquals(1, fromSunId.size());
		assertEquals("d2c35e1", fromSunId.get(0).getScopId());
	}
//	@Test
	public void testHard() {
		MockScopDatabase scop = new MockScopDatabase();
		scop.addClaLineWithHie("d2c35e1	2c35	E:14-139	a.60.8.2	129717	cl=46456,cf=47768,sf=47819,fa=69044,dm=69045,sp=140642,px=129717");
		scop.addDesLine("129717	px	a.60.8.2	d2c35e1	2c35 E:14-139");
	}
}
