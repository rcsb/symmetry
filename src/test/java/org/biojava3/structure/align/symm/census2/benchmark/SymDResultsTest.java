package org.biojava3.structure.align.symm.census2.benchmark;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Alignment;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.protodomain.ResourceList;
import org.biojava3.structure.align.symm.protodomain.ResourceList.NameProvider;
import org.junit.Before;
import org.junit.Test;


public class SymDResultsTest {

	private static final String SYMD_PATH = "src/test/resources/census2/benchmark/symd";

//	@Before
	public void setUp() throws StructureException {
		ResourceList.set(NameProvider.defaultNameProvider(), ResourceList.DEFAULT_PDB_DIR);
	}
	
//	@Test
	public void testRunSymDSimple() throws SymDException {
		final String pdbFile = "src/test/resources/census2/benchmark/1WOP.pdb";
		Result result = SymDResults.runSymD(SYMD_PATH, pdbFile);
		assertEquals("1WOP", result.getScopId());
		final Alignment alignment = result.getAlignment();
		assertEquals(109, (int) alignment.getInitialShift());
		assertEquals(140, (int) alignment.getAlignLength());
		assertEquals(140, (int) alignment.getnNonSelfAligned());
		assertEquals(134.07, (float) alignment.getAlternateTm(), 0.01f);
		assertEquals(0.3683, (float) alignment.getTmpr(), 0.01f);
		assertEquals(10.66, (float) alignment.getzScore(), 0.01f);
	}

//	@Test
	public void testRunSymDMultiple() {
		final String pdbFilesPath = "src/test/resources/census2/benchmark";
		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		final ScopDatabase scop = ScopFactory.getSCOP();
		domains.add(scop.getDomainByScopID("d3ejba1"));
		SymDResults results = SymDResults.runSymD(SYMD_PATH, pdbFilesPath, ResourceList.get().getCache(), domains);
		for (int i = 0; i < domains.size(); i++) {
			final Result result = results.getData().get(i);
			assertEquals(domains.get(i).getScopId(), result.getScopId());
		}
	}
	
//	@Test
	public void testWriteToFile() throws IOException {
		
		File lineByLine = new File("src/test/resources/census2/benchmark/list_for_symd");
		File outputFile = new File("src/test/resources/census2/benchmark/symd_actual_result.xml");
		SymDResults.writeToFile(SYMD_PATH, lineByLine, ResourceList.get().getCache(), outputFile);
		
		String expected = ResourceList.get().openFileAsString("census2/benchmark/symd_expected_result.xml");
	
		// unfortunately, the timestamp will be different
		String[] expectedLines = expected.split("\n");
		BufferedReader br = ResourceList.get().openReader("census2/benchmark/symd_actual_result.xml");
		String line = "";
		int i = 0;
		while ((line = br.readLine()) != null) {
			if (!line.contains("<timestamp>")) {
				assertEquals(expectedLines[i], line);
			}
			i++;
		}
		br.close();
		
		outputFile.delete();
	}
	
}
