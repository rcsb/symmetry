package org.biojava3.structure.align.symm.census2;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.junit.Test;

/**
 * A test for {@AlignmentMapping}.
 * @author dmyersturnbull
 */
public class AlignmentMappingTest {

	@Test
	public void testBuildAfpChain() throws IOException, StructureException {
		Results census = Results.fromXML(new File("src/test/resources/census2/expected1_with_map.xml"));
		Result result = census.getData().get(0);
		AlignmentMapping mapping = result.getAlignmentMapping();
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75);
//		ScopFactory.setScopDatabase(scop);
//		ScopDatabase scop = new BerkeleyScopInstallation();
//		scop.setScopVersion("1.75B");
		ScopDomain domain = scop.getDomainByScopID("d2c35e1");
		System.out.println(domain);
		AtomCache cache = new AtomCache();
		Structure structure = cache.getStructureForDomain(domain);
		Atom[] ca = StructureTools.getAtomCAArray(structure);
		AFPChain afpChain = mapping.buildAfpChain(ca, StructureTools.cloneCAArray(ca));
		assertEquals("Wrong TM-score", 0.24488482, afpChain.getTMScore(), 0.00000001);
	}

}
