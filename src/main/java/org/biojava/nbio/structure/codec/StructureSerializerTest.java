package org.biojava.nbio.structure.codec;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.bio.structure.rcsb.GetRepresentatives;
import org.biojava3.structure.StructureIO;
import org.rcsb.codec.StructureInflator;

public class StructureSerializerTest implements Runnable {
	private AtomCache cache = null;
	private static String RESULT_DIR = "/Users/Peter/Documents/StructureSerializerCache";
	private static String RESULT_DIR_R = "/Users/Peter/Documents/StructureSerializerResults/";
	boolean useFile = false;

	public StructureSerializerTest () {
		initializeCache();
	}

	public static void main(String[] args) {
		new StructureSerializerTest().run();
	}

	@Override
	public void run() {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileWriter(RESULT_DIR_R + timeStamp + "_serializer.csv"));
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}
		out.println("pdbID,io,type,gzipped,compressionLevel,atomCount,fileSize,time");
	
		boolean write = true;
		boolean read = true;
		boolean shuffle = false;
//		boolean useFile = false;
		
		List<String> pdb = new ArrayList<String>(GetRepresentatives.getAll());
//		List<String> pdb = new ArrayList<String>();
//		pdb.add("100D");
		pdb = pdb.subList(0, 100);

		if (shuffle) {
			Collections.shuffle(pdb);
		}
		
		if (write) {
			try {
				writeTests(out, pdb);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//			
		if (read) {
		//	if (shuffle) {
		//		Collections.shuffle(pdb);
		//	}
			try {
				readTests(out, pdb);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	private void readTests(PrintWriter out, List<String> pdb) throws Exception {
		int success = 0;
		int fail = 0;
		
		for (String pdbId: pdb) {
			try {
				int compressionLevel = 1;
				boolean caOnly = false;
				read(out, pdbId, compressionLevel, caOnly);
				success++;
			} catch (IOException e) {
				fail++;
				e.printStackTrace();
			}	
		}

		System.out.println("Structures read success : " + success);
		System.out.println("Structures read failed  : " + fail);
	}

	private void writeTests(PrintWriter out, List<String> pdb) throws Exception {
		int success = 0;
		int fail = 0;
		long totalAtomCount = 0;
//		long totalCAAtomCount = 0;
		long totalTimeGz = 0;
		long totalTimeLz = 0;


		for (String pdbId: pdb) {
			System.out.println("------------- " + pdbId  + "-------------");
			
			try {
				FileParsingParameters params = cache.getFileParsingParams();
				boolean caOnly = false;
				params.setParseCAOnly(caOnly);

				//long t1 = System.nanoTime(); 
				Structure structureAll = StructureIO.getStructure(pdbId);
				//long t2 = System.nanoTime();
				int atomCount = StructureTools.getNrAtoms(structureAll);
				totalAtomCount += atomCount;
				//long time = t2 - t1;
				//			System.out.println("pdb read time: " + (t2-t1) + " atomCount: " + atomCount);
				//		//		out.println("PDB_" + pdbId + "," + "c" + "," + "all" + "," + "TRUE" + "," + "-1" + "," + atomCount + "," + "" + "," + time/1E6);
				//		//		out.flush();
				int compressionLevel = 1;
				//				boolean caOnly = true;
				boolean gzipped = true;
				long start = System.nanoTime();
				write(out, pdbId, structureAll, compressionLevel, caOnly, gzipped);
				long gzTime = System.nanoTime();
				totalTimeGz += (gzTime-start);

				//					read(out, pdbId, compressionLevel, caOnly, gzipped);
				//					totalTimeGz = System.nanoTime() - start;
				success++;
			} catch (Exception e) { // this can be an IOException or StructureException
				fail++;
				e.printStackTrace();
			}
		}

		out.flush();
		
		System.out.println("Structures read success : " + success);
		System.out.println("Structures read failed  : " + fail);
		System.out.println("Total all atom count    : " + totalAtomCount);
//		System.out.println("Total CA atom count     : " + totalCAAtomCount);
		System.out.println("Total gz time           : " + totalTimeGz/1000000);
		System.out.println("Total lz time           : " + totalTimeLz/1000000);
	}

	private void write(PrintWriter out, String pdbId, Structure structure, int compressionMethod, boolean caOnly, boolean gzipped) throws IOException {
		String type = "all";
		if (caOnly) {
			type = "CA";
		}
		String extension = ".hesc";
		String compression = "";
		if (gzipped) {
			extension =".hesc.gz";
			compression = "_gz";
		}
		String io = "w";
		String name = pdbId + extension;
		String fileName = RESULT_DIR + "_" + compressionMethod + "_" + type + compression + "/" + name;

		System.out.println(fileName);
		BioJavaStructureDeflator deflator = new BioJavaStructureDeflator();
		deflator.deflate(structure, fileName, compressionMethod);
		long fileSize = deflator.getFileSizeCompressed();
		long time = deflator.getWriteTime();
		int atomCount = StructureTools.getNrAtoms(structure);
		System.out.println("Write " + fileName + ": " + time/1E6 + " ms");

		out.println("PDB_" + pdbId + "," + io + "," + type + "," + gzipped + "," + compressionMethod + "," + atomCount + "," + fileSize + "," + time/1E6);
	}

	private void read(PrintWriter out, String pdbId, int compressionLevel, boolean caOnly) throws Exception {
		String type = "all";
		if (caOnly) {
			type = "CA";
		}
		String extension = ".hesc";
		String compression = "";
		boolean gzipped = true; // TODO remove gz extension
		if (gzipped) {
			extension =".hesc.gz";
			compression = "_gz";
		}
		String io = "r";
		String name = pdbId + extension;
		String fileName = RESULT_DIR + "_" + compressionLevel + "_" + type + compression + "/" + name;

		BioJavaStructureInflator inflator = new BioJavaStructureInflator();
		StructureInflator def = new StructureInflator(inflator);
		if (useFile) {
			def.read(fileName);
		} else {
		    FileInputStream inputStream = new FileInputStream(fileName);
		    def.read(inputStream);
		}
		int fileSize = (int) def.getFileSizeCompressed();
		long time = def.getReadTime();
		def.close();
		Structure s= inflator.getStructure();
		int atomCount = StructureTools.getNrAtoms(s);
		s = null;
		System.out.println("Read " + fileName + ": " + time/1E6 + " ms, atoms: " + atomCount);

		out.println("PDB_" + pdbId + "," + io + "," + type + "," + gzipped + "," + compressionLevel + "," + atomCount + "," + fileSize + "," + time/1E6);
		out.flush();
	}

	/*
	private void printAtoms(Structure structure) {
		Atom[] atoms = StructureTools.getAllAtomArray(structure);
		for (Atom a: atoms) {
			System.out.println(a.toPDB());
		}
	}
	*/

	private void initializeCache() {
		cache = new AtomCache();
		//		cache.setUseMmCif(true);
		//		System.out.println("cache: " + cache.getPath());
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAtomCaThreshold(Integer.MAX_VALUE);
		params.setAlignSeqRes(true);
		//		params.setParseCAOnly(true);
		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(true);
		//		MmCifBiolAssemblyProvider mmcifProvider = new MmCifBiolAssemblyProvider();
		//		BioUnitDataProviderFactory.setBioUnitDataProvider(mmcifProvider.getClass().getCanonicalName());	
		//		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		StructureIO.setAtomCache(cache);
	}

	//	private static String[] testCase = {"2KU2"};
	//	private static String[] testCase = {"1VU4"}; // nmr core structure seems to be identical
	//	private static String[] testCase = {"2WDK"};
	//	private static String[] testCase = {"2AAZ"}; // large protein, x-ray
	//	private static String[] testCase = {"1HTQ"}; // 10 copies of AU with o=0.1;	
	//	private static String[] testCase = {"1VU4"};
	//	private static String[] testCase = {"2M8R"};
	//	private static String[] testCase = {"4HHB"};	
	//	private static String[] testCase = {"1STP"};
	//	private static String[] testCase = {"1AO2"};
	//	private static String[] testCase = {"1BPV"};
	//	private static String[] testCase = {"1CDG"};
	//	private static String[] testCase = {"3NHD"};
	//	private static String[] testCase = {"1STP"};
	//	private static String[] testCase = {"3j3q"};
	//	private static String[] testCase = {"2K8M","1HTQ","1STP","2KU2","1VU4"};
	// static String[] testCase = {"1VU4","2AAZ","4HHB","2KU2","1BPV","2WDK"};
//	static String[] testCase = {"1VU4","2AAZ","4HHB","2KU2","1BPV","2WDK","1HTQ"};
//	static String[] testCase = {"1STP","3NHD","3LOD","2GJI","2K8M","1VU4","4HHB","2KU2","1BPV","2WDK","1HTQ"};
	//	private static String[] testCase = {"1E3M"};

	//	Map<Integer, List<BiologicalAssemblyTransformation>> bas = structure.getPDBHeader().getBioUnitTranformationMap();
	//					Date date = structure.getPDBHeader().getDepDate();
	//					System.out.println(date.toString());
	//					structure.getCrystallographicInfo();
	//					structure.getName();
	//					structure.getPDBHeader().getTitle();
	//					structure.getPDBHeader().getResolution();
	//					structure.getPDBHeader().getTechnique();
}
