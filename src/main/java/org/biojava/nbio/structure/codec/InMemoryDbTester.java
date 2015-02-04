package org.biojava.nbio.structure.codec;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.rcsb.GetRepresentatives;
import org.biojava.nbio.structure.StructureIO;


public class InMemoryDbTester implements Runnable {
	private static String RESULT_DIR_R = "C:/Users/Peter/Documents/StructureSerializerResults/";
	
	public static void main(String[] args) {
		new InMemoryDbTester().run();
	}

	public void run() {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		PrintWriter out = null;

		try {
			out = new PrintWriter(new FileWriter(RESULT_DIR_R + timeStamp + "_database.csv"));
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}

		out.println("pdbID,atomCount,fileSize,fileSizeCompressed,compressionLevel,pdbTime,writeCompressed,writeUncompressed," + 
				"readCompressed,readUncompressed");


		List<String> pdb = new ArrayList<String>(GetRepresentatives.getAll());
		
		pdb = pdb.subList(0, 10000);

		long pdbTime = 0;
		long fileSize = 0;
		long fileSizeCompressed = 0;
		long writeCompressed = 0;
		long writeUncompressed = 0;
		long readCompressed = 0;
		long readUncompressed = 0;
		int count = 0;
		long totalAtomCount = 0;
		int useCase = 1;
		int compressionLevel = 2;
//		Collections.shuffle(pdb);
		
		List<String> testCases = Arrays.asList(testCase);
		
		
		long start = System.nanoTime();
		RcsbPdbInMemoryDatabase db = new RcsbPdbInMemoryDatabase();
		try {
	//		db.addPdbIds(testCases);
			db.addPdbIds(pdb);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long end = System.nanoTime();
		System.out.println("Database creation time: " + (end - start)/1E6 + " ms");
		System.out.println("Database size: " + db.getSize());
     
//		count = 0;
//		long totalReadTime = 0;
//		long t1 = System.nanoTime();
//		for (int i = 0; i < pdb.size()-1; i++) {
//			System.out.println("------------- " + pdb.get(i)  + "-------------");
//			BioJavaPdbDeflator deflator = new BioJavaPdbDeflator();
//			try {
//				db.getData(pdb.get(i), deflator);
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//			Structure s1 = deflator.getStructure();
//			for (int j = i+1; j < pdb.size(); j++) {
//				BioJavaPdbDeflator deflator2 = new BioJavaPdbDeflator();
//				try {
//					db.getData(pdb.get(j), deflator2);
//				} catch (IOException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
//				Structure s2 = deflator.getStructure();
//				s2 = null;
//			}
//			s1 = null;
//		}
//		totalReadTime += System.nanoTime() - t1;

		Collections.shuffle(pdb);
		count = 0;
		long totalReadTime = 0;

	//	for (String pdbId: testCases) {

		for (String pdbId: pdb) {

			System.out.println("------------- " + pdbId  + "-------------");
			
			long t1 = System.nanoTime();
//			BioJavaPdbDeflator deflator = new BioJavaPdbDeflator();
//			try {
//				db.getData(pdbId, deflator);
//				count++;
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			} catch (Exception e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
			totalReadTime += System.nanoTime() - t1;
			//Structure structure = deflator.getStructure();
			//int atomCount = StructureTools.getNrAtoms(structure);
			//System.out.println("atomCount: " + atomCount);
 


			//			out.println("PDB_" + pdbId + "," + atomCount + "," + fileSize + "," + fileSizeCompressed + "," + compressionLevel +  "," + pdbTime/1E6 + ","  
			//					+ writeCompressed/1E6 + "," + writeUncompressed/1E6 + "," + readCompressed/1E6 + "," + readUncompressed/1E6);
			//			out.flush();
		}
		
		System.out.println("Database creation time: " + (end - start)/1E9 + " s");
		System.out.println("Total atom count: " + totalAtomCount);
		System.out.println("Structure read: " + count);
		System.out.println("Total read time: " + totalReadTime/1E9 + " s");
	}
	
	private void printAtoms(Structure structure) {
		Atom[] atoms = StructureTools.getAllAtomArray(structure);
		for (Atom a: atoms) {
			System.out.println(a.toPDB());
		}
	}

	
//	private static String[] testCase = {"2KU2"}; // nmr core structure seems to be identical
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
//	private static String[] testCase = {"13PK"};
 static String[] testCase = {"1VU4","2AAZ","4HHB","2KU2","1BPV","2WDK"};

}
