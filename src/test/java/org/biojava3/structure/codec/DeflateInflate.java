package org.biojava3.structure.codec;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import org.biojava.bio.structure.Atom;
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
		String[] pdbIds = {"1ZMP"};
//		String[] pdbIds = {"1STP","4HHB","1OHR"};

		for (String pdbId: pdbIds) {
			Structure original = getStructure(pdbId);
			String fileName = deflate(original, pdbId);
			Structure copy = inflate(fileName);

			int expectedCount =  StructureTools.getNrAtoms(original);
			int actualCount = StructureTools.getNrAtoms(copy);	
			assertEquals(expectedCount, actualCount);

			Atom[] expectedAtoms = StructureTools.getAllAtomArray(original);
			Atom[] actualAtoms = StructureTools.getAllAtomArray(copy);

			for (int i = 0; i < expectedAtoms.length; i++) {
	//			System.out.println(expectedAtoms[i].toPDB());
	//			System.out.println(actualAtoms[i].toPDB());
				assertEquals(maskSerialNumber(expectedAtoms[i].toPDB()), maskSerialNumber(actualAtoms[i].toPDB()));
			}
		}
	}
	
	public static String deflate(Structure structure, String pdbId) throws IOException {
		File temp = File.createTempFile(pdbId, CodecConstants.FileExtension);
		String fileName = temp.getName();
		int compressionLevel = 1;
		
		BioJavaStructureDeflator deflator = new BioJavaStructureDeflator();
		deflator.deflate(structure, fileName, compressionLevel);
		
		return fileName;
	}
	
	public static Structure inflate(String fileName) throws Exception {
		BioJavaStructureInflator inflator = new BioJavaStructureInflator();
		StructureInflator def = new StructureInflator(inflator);
	    FileInputStream inputStream = new FileInputStream(fileName);
		def.read(inputStream);
		return inflator.getStructure();
	}
	
	private static Structure getStructure(String pdbId) throws IOException, StructureException {
		initializeCache();
		return StructureIO.getStructure(pdbId);
	}
	
	private String maskSerialNumber(String atomRecord) {
		return atomRecord.substring(0,  6) + "     " + atomRecord.substring(11, atomRecord.length());
	}
	
	private static void initializeCache() {
		AtomCache cache = new AtomCache();
		//		System.out.println("cache: " + cache.getPath());
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
