package demo;

import static org.junit.Assert.assertEquals;

import java.io.FileNotFoundException;
import java.net.URL;
import java.util.List;

import org.junit.Test;

import demo.CeSymmMain;

public class CeSymmMainTest {

	@Test
	public void testParseInput() throws FileNotFoundException {
		URL url = CeSymmMainTest.class.getResource("/cesymmtest.txt");
		String filename = url.getFile();
		List<String> names = CeSymmMain.parseInputStructures(filename);
		
		assertEquals(6,names.size());
		int i = 0;
		assertEquals("d1ijqa1",names.get(i++));
		assertEquals("1G6S",names.get(i++));
		assertEquals("1MER.A",names.get(i++));
		assertEquals("d1h70a_",names.get(i++));
		assertEquals("2YMS_A:,C:,B:,D:",names.get(i++));
		assertEquals("1HIV",names.get(i++));
	}

}
