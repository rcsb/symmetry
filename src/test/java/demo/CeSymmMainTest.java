package demo;

import static org.junit.Assert.assertEquals;

import java.io.FileNotFoundException;
import java.net.URL;
import java.util.List;

import org.junit.Test;

public class CeSymmMainTest {

	@Test
	public void testParseInput() throws FileNotFoundException {
		URL url = CeSymmMainTest.class.getResource("/inputtest.txt");
		String filename = url.getFile();
		List<String> names = CeSymmMain.parseInputStructures(filename);
		
		assertEquals(6,names.size());
		int i = 0;
		assertEquals("line1",names.get(i++));
		assertEquals("line2.1",names.get(i++));
		assertEquals("2.2",names.get(i++));
		assertEquals("twothree",names.get(i++));
		assertEquals("three",names.get(i++));
		assertEquals("four",names.get(i++));
	}

}
