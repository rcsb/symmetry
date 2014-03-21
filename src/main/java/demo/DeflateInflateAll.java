package demo;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.codec.BioJavaStructureDeflator;
import org.biojava3.structure.codec.BioJavaStructureInflator;
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.rcsb.codec.CodecConstants;
import org.rcsb.codec.StructureInflator;

public class DeflateInflateAll {


	public void test() throws Exception {
		List<String> pdbIds = new ArrayList<String>(GetRepresentatives.getAll());
		pdbIds.remove("11GS"); // PRO A-2 listed twice in seqres groups (issue with alignment)
		pdbIds.remove("12GS"); // PRO A-2 listed twice in seqres groups (issue with alignment)
		pdbIds.remove("13GS"); // PRO A-2 listed twice in seqres groups (issue with alignment)
		pdbIds.remove("14GS"); // PRO A-2 listed twice in seqres groups (issue with alignment)
		pdbIds.remove("16GS"); // PRO A-2 listed twice in seqres groups (issue with alignment)
		pdbIds.remove("17GS"); // PRO A-2 listed twice in seqres groups (issue with alignment)
		pdbIds.remove("18GS"); // PRO A-2 listed twice in seqres groups (issue with alignment)
		pdbIds.remove("136D"); // issue with irregular numbering of inserted residues, or missing seq. res. residues in model 1??
		pdbIds.remove("177D"); // multiple model issue
		pdbIds.remove("176D"); // B: GAGUUC: UUC added twice, once without atoms, once with atoms? Problem handling nucleotides??
		pdbIds.remove("148L"); // non-std. amino acids (D, isopeptide?) are part of chain
		pdbIds.remove("1A07"); // hetatm in chain C ??
		pdbIds.remove("1A08"); // hetatm in chain C ??
		pdbIds.remove("1A09"); // hetatm in chain C ??
		pdbIds.remove("1A1A"); // hetatm in chain C ??
		pdbIds.remove("1A1B"); // hetatm in chain C ??
		pdbIds.remove("1A1C"); // hetatm in chain C ??
		pdbIds.remove("1A1E"); // hetatm in chain C ??
		pdbIds.remove("1E3M"); // has missing chain id in link records -> StructureException; issue reported to RU
		pdbIds.remove("1GVX"); // invalid link record
		pdbIds.remove("1OAO"); // ..
		pdbIds.remove("1QJH"); // ..
		pdbIds.remove("1QJI"); // ..
		pdbIds.remove("1GUG"); // issue with alt loc in link record
	    pdbIds.remove("1KO5"); // ..
	    pdbIds.remove("1LLB"); // ..
//	    pdbIds.remove("1NPQ"); // file has  -0.000 vs. 0.000
//	    pdbIds.remove("1PUL"); // file has  -0.000 vs. 0.000
//	    pdbIds.remove("1Q8K"); // file has  -0.000 vs. 0.000
	    pdbIds.remove("1OAX"); // could not find chain "I" ")
	    pdbIds.remove("1OAY"); // could not find chain "I" ")
	    pdbIds.remove("1X26"); // HETATM       C1 NNAZ A vs. HETATM       C1  NAZ A  25
	    // compression/decompression errors:
	    pdbIds.remove("2V93"); // 	for (int k = 0; k < atomCount; k++) {String atomName = info[index++]; index of of bounds?
	    // this is caused by getName() in AtomImpl: a = alt.getAtom(name);
		// dirty hack
		// we are adding this group to the main one...
	    
	    pdbIds.remove("2Y8V"); // ATOM         CB  GLU C  52      30.579  68.881  95.567  1.00 21.01           C vs ATOM         CB  GLU C  52      30.579  68.881 100.366  1.00 21.01           C
        pdbIds.remove("3K1Q"); // ATOM         CG  HIS Y  37 vs. ATOM       0  CG  HIS Y  37

		boolean skip = false;
		String startId = "3K1Q";

//		pdbIds = Arrays.asList("1HRH");
		
		for (String pdbId: pdbIds) {
			if (pdbId.equals(startId)) {
				skip = false;
			}
			if (skip) continue;
			System.out.println(pdbId);
			System.out.println("---------------" + pdbId + "----------------");
			Structure original = null;
			try {
				original = getStructure(pdbId);
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			} 
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
		File temp = File.createTempFile(pdbId, CodecConstants.CODEC_FILE_EXTENSION);
		String fileName = temp.getName();
		int compressionMethod = 1;
		
		BioJavaStructureDeflator deflator = new BioJavaStructureDeflator();
		deflator.deflate(structure, fileName, compressionMethod);
		
		return fileName;
	}
	
	public static Structure inflate(String fileName) throws Exception {
		BioJavaStructureInflator inflator = new BioJavaStructureInflator();
		StructureInflator def = new StructureInflator(inflator);
	    FileInputStream inputStream = new FileInputStream(fileName);
		def.read(inputStream);
		inputStream.close();
		return inflator.getStructure();
	}
	
	private static Structure getStructure(String pdbId) throws IOException, StructureException {
		initializeCache();
		return StructureIO.getStructure(pdbId);
	}
	
	private String maskSerialNumber(String atomRecord) {
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
		index = sb.lastIndexOf("-0.00"); // negative zero of b-factors TODO this could also match a coordinate
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

	private void assertEquals(String maskSerialNumber, String maskSerialNumber2) {
		if (! maskSerialNumber.equals(maskSerialNumber2))
			throw new RuntimeException(maskSerialNumber + " != " + maskSerialNumber2);
		
	}

	private void assertEquals(int expectedCount, int actualCount) {
		if (expectedCount != actualCount)
			throw new RuntimeException(expectedCount + " != " + actualCount);
		
	}

}
