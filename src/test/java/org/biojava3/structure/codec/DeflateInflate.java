package org.biojava3.structure.codec;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava3.structure.StructureIO;
import org.junit.Test;
import org.rcsb.codec.CodecConstants;
import org.rcsb.codec.StructureInflator;

public class DeflateInflate {

	@Test
	public void test() throws Exception {
//		String[] pdbIds = {"1STP"};
		String[] pdbIds = {"18GS"};
//		String[] pdbIds = {"148L"};
//		String[] pdbIds = {"11GS"};
//		String[] pdbIds = {"1STP","4HHB","1OHR"};

		for (String pdbId: pdbIds) {
			Structure original = getStructure(pdbId);
			
			printFirstTwoResidues(original);
			
			File file = deflate(original, pdbId);
			Structure copy = null;
			try {
				copy = inflate(file);
			} catch (Exception e){
				file.delete();
				e.printStackTrace();

				fail(e.getMessage());

			}
			file.delete();
			// it is just a tmp file, clean up..
			assertNotNull(copy);
			
			
			int expectedCount =  StructureTools.getNrAtoms(original);
			int actualCount = StructureTools.getNrAtoms(copy);	
			assertEquals("Original structure and copy don't have the same number of atoms!" ,expectedCount, actualCount);

			Atom[] expectedAtoms = StructureTools.getAllAtomArray(original);
			Atom[] actualAtoms = StructureTools.getAllAtomArray(copy);

			for (int i = 0; i < expectedAtoms.length; i++) {
	//			System.out.println(expectedAtoms[i].toPDB());
	//			System.out.println(actualAtoms[i].toPDB());
				assertEquals(maskSerialNumber(expectedAtoms[i].toPDB()), maskSerialNumber(actualAtoms[i].toPDB()));
			}
		}
	}
	
	private void printFirstTwoResidues(Structure original) {
		
		Chain c = original.getChain(0);
		
		Group first  = c.getSeqResGroup(0);
		Group second = c.getSeqResGroup(1);
		Group third  = c.getSeqResGroup(2);
		
		System.out.println(first);
		System.out.println(second);
		System.out.println(third);
		
	}

	public static File deflate(Structure structure, String pdbId) throws IOException {
		File temp = File.createTempFile(pdbId, CodecConstants.CODEC_FILE_EXTENSION);

		String fileName = temp.getPath();
		System.out.println(fileName);
		int compressionMethod = 1;
		
		BioJavaStructureDeflator deflator = new BioJavaStructureDeflator();
		deflator.deflate(structure, fileName, compressionMethod);
		System.out.println("Compressed file size: " + deflator.getFileSizeCompressed());
		
		return temp;
	}
	
	public static Structure inflate(File file) throws Exception {
		BioJavaStructureInflator inflator = new BioJavaStructureInflator();
		StructureInflator def = new StructureInflator(inflator);
	    FileInputStream inputStream = new FileInputStream(file);
		def.read(inputStream);
		return inflator.getStructure();
	}
	
	private static Structure getStructure(String pdbId) throws IOException, StructureException {
		initializeCache();
		return StructureIO.getStructure(pdbId);
	}
	
	private String maskSerialNumber(String atomRecord) {
//		System.out.println(atomRecord);
		atomRecord = replaceMinusZero(atomRecord);
		return atomRecord.substring(0,  6) + "     " + atomRecord.substring(11, atomRecord.length());
	}
	
	private String replaceMinusZero(String atomRecord) {
		StringBuffer sb = new StringBuffer(atomRecord);
		int index = sb.indexOf("-0.000"); // negative zero of coordinates
		while (index > 0) {
			sb.setCharAt(index, ' ');
			index = sb.indexOf("-0.000");
		}
		index = sb.lastIndexOf("-0.00");
		if (index >= 60) {
			sb.setCharAt(index, ' ');
		}
		return sb.toString();
	}
	
	private static void initializeCache() {
		AtomCache cache = new AtomCache();

		System.out.println("cache: " + cache.getPath());
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAtomCaThreshold(Integer.MAX_VALUE);
		params.setAlignSeqRes(true);
		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(true);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		StructureIO.setAtomCache(cache);
	}


}
