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
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.junit.Test;

/**
 * A test for {@AlignmentMapping}.
 * @author dmyersturnbull
 */
public class AlignmentMappingTest {

	public static void main(String[] args) throws IOException {
		AtomCache cache = new AtomCache();
		ScopFactory.setScopDatabase(ScopFactory.LATEST_VERSION);
		CensusJob job = CensusJob.setUpJob("d1ux8a_", 0, Census.AlgorithmGiver.getDefault(), Census.getDefaultSignificance(), cache, ScopFactory.getSCOP());
		job.setRecordAlignmentMapping(true);
		Result result = job.call();
		Results results = new Results();
		results.add(result);
		System.out.println(results.toXML());
	}
	
	@Test
	public void testBuildAfpChain() throws IOException, StructureException {
		Results census = Results.fromXML(new File("src/test/resources/census2/expected1_with_map.xml"));
		Result result = census.getData().get(0);
		AlignmentMapping mapping = result.getAlignmentMapping();
		ScopDomain domain = ScopFactory.getSCOP().getDomainByScopID("d1ux8a_");
		AtomCache cache = new AtomCache();
		Structure structure = cache.getStructure(domain.getScopId());
		Atom[] ca1 = StructureTools.getAtomCAArray(structure);
		Atom[] ca2 = StructureTools.getAtomCAArray(structure);
		AFPChain afpChain = mapping.buildAfpChain(ca1, ca2);
		assertEquals("Wrong TM-score", 0.509579, afpChain.getTMScore(), 0.00000001);
	}

}
