package script;



import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.bio.structure.quaternary.io.BioUnitDataProviderFactory;
import org.biojava.bio.structure.quaternary.io.MmCifBiolAssemblyProvider;
import org.biojava.bio.structure.quaternary.io.PDBBioUnitDataProvider;
import org.biojava3.structure.StructureIO;
import org.rcsb.fatcat.server.PdbChainKey;
import org.rcsb.fatcat.server.dao.SequenceClusterDAO;

public class CompareBioAssemblyProviders {
	public static void main(String[] args){
		
		StructureIO.setPdbPath("/Volumes/Macintosh HD2/PDB/");
		
		CompareBioAssemblyProviders me = new CompareBioAssemblyProviders();

		me.run();

		//me.testID("1A02", 1);
	}

	public void run(){

		

		SequenceClusterDAO dao = new SequenceClusterDAO();

		SortedSet<PdbChainKey> representatives = dao.getClusterEntities(90);

		System.out.println("got " + representatives.size() + " representatives...");

		List<String> testedIds = new ArrayList<String>();
		int count = 0;
		for ( PdbChainKey key : representatives){

			String pdbID = key.getPdbId();

			if ( testedIds.contains(pdbID))
				continue;

			count++;

			System.out.println("#" + count + " " + pdbID);

			testedIds.add(pdbID);
			if ( StructureIO.hasBiologicalAssembly(pdbID))
				testID(pdbID,1);

		}
	}


	private void testID(String pdbId, int bioMolecule){


		try {
			// get bio assembly from PDB file
			PDBBioUnitDataProvider pdbProvider = new PDBBioUnitDataProvider();
			BioUnitDataProviderFactory.setBioUnitDataProvider(pdbProvider.getClass().getCanonicalName());
			Structure pdbS = StructureIO.getBiologicalAssembly(pdbId, bioMolecule);

			// get bio assembly from mmcif file
			MmCifBiolAssemblyProvider mmcifProvider = new MmCifBiolAssemblyProvider();
			BioUnitDataProviderFactory.setBioUnitDataProvider(mmcifProvider.getClass().getCanonicalName());			
			Structure mmcifS = StructureIO.getBiologicalAssembly(pdbId, bioMolecule);

			BioUnitDataProviderFactory.setBioUnitDataProvider(BioUnitDataProviderFactory.DEFAULT_PROVIDER_CLASSNAME);



			PDBHeader pHeader = pdbS.getPDBHeader();
			PDBHeader mHeader = mmcifS.getPDBHeader();
			//PDBHeader fHeader = flatFileS.getPDBHeader();

			assertTrue("not correct nr of bioassemblies " + pHeader.getNrBioAssemblies() + " " , pHeader.getNrBioAssemblies() >= bioMolecule);
			assertTrue("not correct nr of bioassemblies " + mHeader.getNrBioAssemblies() + " " , mHeader.getNrBioAssemblies() >= bioMolecule);
			//assertTrue("not correct nr of bioassemblies " + fHeader.getNrBioAssemblies() + " " , fHeader.getNrBioAssemblies() >= bioMolecule);

			// mmcif files contain sometimes partial virus assemblies, so they can contain more info than pdb
			assertTrue(pHeader.getNrBioAssemblies() <= mHeader.getNrBioAssemblies());


			Map<Integer, List<BiologicalAssemblyTransformation>> pMap = pHeader.getBioUnitTranformationMap();
			Map<Integer, List<BiologicalAssemblyTransformation>> mMap = mHeader.getBioUnitTranformationMap();

			//System.out.println("PDB: " + pMap);

			//System.out.println("Mmcif: " + mMap);

			assertTrue(pMap.keySet().size()<= mMap.keySet().size());

			//			for ( Integer k : pMap.keySet()) {
			//				assertTrue(mMap.containsKey(k));
			//				
			//				List<ModelTransformationMatrix> pL = pMap.get(k);
			//				
			//				// mmcif list can be longer due to the use of internal chain IDs
			//				List<ModelTransformationMatrix> mL = mMap.get(k);
			//				
			//				//assertEquals(pL.size(), mL.size());
			//				
			//				
			//				for (ModelTransformationMatrix m1 : pL){
			//					
			//					boolean found = false;
			//					for ( ModelTransformationMatrix m2 : mL){
			//						
			//						if  (! m1.getNdbChainId().equals(m2.getNdbChainId()))
			//								continue;
			//						if ( ! m1.getMatrix().toString().equals(m2.getMatrix().toString()))
			//								continue;
			//						if ( ! equalVectors(m1.getVector(),m2.getVector()))
			//							continue;
			//						
			//						found = true;
			//						
			//					}
			//					
			//					if ( ! found ){
			//						System.err.println("did not find matching matrix " + m1);
			//						System.err.println(mL);
			//					}
			//					assertTrue(found);
			//					
			//				}
			//			}


			assertEquals("Not the same number of chains!" , pdbS.size(),mmcifS.size());

			Atom[] pdbA = StructureTools.getAllAtomArray(pdbS);

			Atom[] mmcifA = StructureTools.getAllAtomArray(mmcifS);

			assertEquals(pdbA.length, mmcifA.length);

	
			// the atom pos can change...
			assertEquals(pdbA[0].toPDB(), mmcifA[0].toPDB());


			// compare with flat file version:
			AtomCache cache = new AtomCache();
			FileParsingParameters params = cache.getFileParsingParams();
			params.setAlignSeqRes(true);
			params.setParseCAOnly(false);

			Structure flatFileS = cache.getBiologicalAssembly(pdbId, bioMolecule, false);

			Atom[] fileA = StructureTools.getAllAtomArray(flatFileS);

			assertEquals(pdbA.length, fileA.length);

		} catch (Exception e){
			e.printStackTrace();
			throw new RuntimeException(e.getMessage());
		}





	}
	private void assertTrue(String string, boolean b) {
		if ( ! b)
			throw new RuntimeException(string);

	}

	private void assertEquals(String string, int size, int size2) {
		if ( size != size2)
			throw new RuntimeException(string);

	}

	private void assertTrue(boolean found) {
		if ( ! found)
			throw new RuntimeException("not true");
	}

	private void assertEquals(String pdb, String pdb2) {
		if ( ! pdb.substring(20).equals(pdb2.substring(20))){

			for (int i = 0; i < pdb.length(); i++) {

				if ( i > pdb2.length())
					break;
				
				if (pdb.charAt(i) != pdb2.charAt(i)) {
					
					System.err.println("first mismatch at pos: " + i + " " + pdb.charAt(i));
					break;
				}
			}

			System.err.println(pdb);
			System.err.println(pdb2);
			System.err.println("Strings don't match");
		}
	}

	private void assertEquals(int length, int length2) {
		if ( length != length2)
			throw new RuntimeException("length don't match");

	}

	private boolean equalVectors(double[] vector, double[] vector2) {

		String s1 = String.format("%.5f %.5f %.5f", vector[0], vector[1], vector[2]);
		String s2 = String.format("%.5f %.5f %.5f", vector2[0], vector2[1], vector2[2]);
		//System.out.println(s1 + " " + s2);
		return s1.equals(s2);

	}
}
