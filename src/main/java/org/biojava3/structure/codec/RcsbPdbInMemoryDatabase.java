/**
 * 
 */
package org.biojava3.structure.codec;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Peter
 *
 */
public class RcsbPdbInMemoryDatabase {
//	private static String RESULT_DIR_C = "C:/Users/Peter/Documents/StructureSerializer_8_all_gz/";
	private static String RESULT_DIR_C = "C:/Users/Peter/Documents/StructureSerializer_8_all_";
	private String format = null;
	private Map<String,byte[]> pdbs = new HashMap<String,byte[]>();
	
	public RcsbPdbInMemoryDatabase() {
		this.format = "gz";
	}
	public RcsbPdbInMemoryDatabase(String format) {
		this.format = format;
	}

	public void addPdbIds(List<String> pdbIds) throws IOException {
		for (String pdbId: pdbIds) {
			System.out.println("Adding: " + pdbId);
			String fileNameCompressed = RESULT_DIR_C + format + "/" + pdbId + ".ser." + format;			
			File file = new File(fileNameCompressed);
			if (file.exists()) {
				byte[] data = new byte[(int) file.length()];
				DataInputStream dis = new DataInputStream(new FileInputStream(file));
				dis.readFully(data);
				dis.close();
				pdbs.put(pdbId, data);
				data = null;
			}
			// Java 7
			// byte[] bytes = Files.readAllBytes(Paths.get(filename));
		}
	}
	
	public int getSize() {
		return pdbs.size();
	}
//	public void getData(String pdbId, RcsbPdbStructureInterface deflator) throws Exception {
//		RcsbPdbDeflator def = new RcsbPdbDeflator(deflator);
//		byte[] data = pdbs.get(pdbId);
//		def.read(data, format);
//		def.close();
//	}
}
