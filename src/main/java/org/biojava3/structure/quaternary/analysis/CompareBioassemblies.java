package org.biojava3.structure.quaternary.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.bio.structure.quaternary.io.BioUnitDataProviderFactory;
import org.biojava.bio.structure.quaternary.io.MmCifBiolAssemblyProvider;
import org.biojava3.structure.StructureIO;
import org.biojava3.structure.dbscan.GetRepresentatives;

public class CompareBioassemblies {
	private AtomCache cache = null;
	private static String RESULT_DIR = "C:/Users/Peter/Documents/QuatStructureComparison/";


	public CompareBioassemblies() {
		initializeCache();
	}

	public static void main(String[] args) {
		new CompareBioassemblies().run();
	}

	public void run() {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		StructureIO.setAtomCache(cache);

		PrintWriter out = null;
		PrintWriter error = null;

		try {
			out = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_report.txt"));
			error = new PrintWriter(new FileWriter(RESULT_DIR + timeStamp + "_error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}

		int success = 0;
		int failure = 0;
		int structureNull = 0;

		Set<String> pdbAll = GetRepresentatives.getAll();

		// set skip to true to restart calculation with a specified PDB ID
		boolean skip = false;
		String restartId = "1A6S";

		for (String pdbId: pdbAll) {
//		for (String pdbId: testCases) {
			if (skip && pdbId.equals(restartId)) {
				skip = false;
			} 
			if (skip) {
				continue;
			}

			if (pdbId.equals("1M4X")) {
				continue;
			}

			int bioAssemblyCount = StructureIO.getNrBiologicalAssemblies(pdbId);
			System.out.println("BioUnitDataProvider: " + BioUnitDataProviderFactory.getBioUnitDataProvider().getClass().getName());
			System.out.println(pdbId + " Bioassemblies: " + bioAssemblyCount);


			
			for (int i = 0; i <= bioAssemblyCount; i++) {	
				System.out.println("------------- " + pdbId  + "[" + i + "] -------------");
		    	System.out.println(ChemCompGroupFactory.getChemCompProvider().getClass().getName());
				long t1 = System.nanoTime();
		    	Structure s1 = createBioAssembly(error, pdbId, i);
		    	System.out.println(ChemCompGroupFactory.getChemCompProvider().getClass().getName());
		    	long t2 = System.nanoTime();   	
		    	Structure s2 = readBioAssemblyFile(error, pdbId, i);
		    	System.out.println(ChemCompGroupFactory.getChemCompProvider().getClass().getName());
		    	long t3 = System.nanoTime();
		    	
		    	System.out.println(pdbId + " On the fly: " + (t2-t1)/1000000 + " From file: " + (t3-t2)/1000000);
		    	
		    	if (s1 == null || s2 == null) {
		    		System.out.println("Structure is null");
		    		structureNull++;
		    		continue;
		    	}
		    	
                String result = compareStructures(s1, s2, pdbId, i);
                if (!result.isEmpty()) {
                	 String prefix = (i > 0) ? pdbId + "[" + i + "]: " : pdbId + ": ";
                     out.println(prefix + result);
                     out.flush();
                     failure++;
                } else {
                	success++;
                }
			}	
    
		}

		System.out.println("PDBs succeeded: " + success);
		System.out.println("PDBs failed   : " + failure);
		System.out.println("PDBs  null    : " + structureNull);
		System.out.println("Total bioassemblies : " + (success + failure+ structureNull));
		System.out.println("Total structure: " + pdbAll.size());

		out.flush();
		out.close();
        error.flush();
		error.close();
	}
	
	private String compareStructures(Structure s1, Structure s2, String pdbId, int bioAssemblyId) {
		Atom[] ca1 = StructureTools.getAtomCAArray(s1);	
		Atom[] ca2 = StructureTools.getAtomCAArray(s2);
		
//		System.out.println("ca1: " + Arrays.toString(ca1));
//		System.out.println("ca2: " + Arrays.toString(ca2));
		if (ca1.length != ca2.length && ! s1.isNmr()) {
			return "Inconsistent number of Calpha atoms: " + ca1.length + " - " + ca2.length;
		}
		
		// skip inconsistent NMR models for now
		if (s1.nrModels() < s2.nrModels() && s1.isNmr()) {
			return "";
		}
		
		String sb1 = getModelChainString(s1);
		String sb2 = getModelChainString(s2);
		System.out.println("Chain sequence: " + sb1 + " - " + sb2);
		if (! sb1.equals(sb2)) {

			return "Inconsistent chain sequence: " + sb1 + " - " + sb2;
		}
		
		
		for (Atom a1: ca1) {
			boolean match = false;
			for (Atom a2: ca2) {
				if (compareDoubleArray(a1.getCoords(), a2.getCoords())) {
//				if (a1.getGroup().getPDBName().equals(a2.getGroup().getPDBName()) &&
//						a1.getElement() == a2.getElement() && 
//						a1.getName().equals(a2.getName()) && 	
//						compareDoubleArray(a1.getCoords(), a2.getCoords())) {
					match = true;
					break;
				}
				
			}
			if (! match) {
				return "Cannot match atom: " + a1;
			}
		}
		
        
        // System.out.println("Models: " + s1.nrModels());
		if (s1.nrModels() != s2.nrModels() && ! s1.isNmr()) {
	      	return "Inconsistent number of models: " + s1.nrModels() + " - " + s2.nrModels();
	    }
        
        for (int i = 0; i < s1.nrModels(); i++) {
        	List<Chain> chains1 = s1.getChains(i);
        	List<Chain> chains2 = s2.getChains(i);
        	if (chains1.size() != chains2.size()) {
        		return "Inconsistent number of chains: " + chains1.size() + " - " + chains2.size();
        	}
        	// System.out.println("Number of chains: " + chains1.size());
//        	for (int j = 0; j < chains1.size(); j++) {
//        		Chain c1 = chains1.get(j);
//        		Chain c2 = chains2.get(j);
//        		if (!c1.getChainID().equals(c2.getChainID())) {
//        			return "Inconsistent chain ids: " + c1.getChainID() + " - " + c2.getChainID();
//        		}
//        		List<Group> groups1 = c1.getAtomGroups();
//        		List<Group> groups2 = c2.getAtomGroups();
//        		if (groups1.size() != groups2.size()) {
//        			return "Inconsistent number of atom groups: " + groups1.size() + " - " + groups2.size();
//        		}
//        		for (int k = 0; k < groups1.size(); k++) {
//        			Group g1 = groups1.get(k);
//        			Group g2 = groups2.get(k);
//        			if (!g1.getPDBName().equals(g2.getPDBName())) {
//        				System.out.println("Group: " + g1);
//   //     				return "Inconsistent PDB names: " + g1.getPDBName() + " - " + g2.getPDBName();
//        			}
//        			List<Atom> atoms1 = g1.getAtoms();
//        			List<Atom> atoms2 = g2.getAtoms();
//        			if (atoms1.size() !=atoms2.size()) {
//        				return "Inconsistent number of atoms: " + atoms1.size() + " - " + atoms2.size();
//        			}
//        			for (int m = 0; m < atoms1.size(); m++) {
//        				Atom a1 = atoms1.get(m);
//        				Atom a2 = atoms2.get(m);
//        				if (a1.getElement() != a2.getElement()) {
//        					return "Inconsistent elements: " + a1.getElement() + " - " + a2.getElement();
//        				}
//        				if (!compareDoubleArray(a1.getCoords(), a2.getCoords())) {
//        					return "Inconsistent coordinates: " + a1 + " - " + a2;
//        				}
//        				// System.out.println(a1);
//        				// System.out.println(a2);
//        			}
//        		}
//        	}
        }
        return "";
	}
	
	private String getModelChainString(Structure s) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < s.nrModels(); i++) {
			sb.append(i);
		    for (Chain c: s.getChains(i)) {
		    	sb.append(c.getChainID());
		    }
		}
		return sb.toString();
	}

	private boolean compareDoubleArray(double[] a1, double[] a2) {
		for (int i = 0; i < a1.length; i++) {
			if (Math.abs(a1[i]-a2[i]) > 0.001) {
				return false;
			}
		}
		return true;
	}


	private void initializeCache() {
		cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(true);
		params.setAtomCaThreshold(Integer.MAX_VALUE);
		params.setLoadChemCompInfo(true);
		params.setMaxAtoms(Integer.MAX_VALUE);
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
		MmCifBiolAssemblyProvider mmcifProvider = new MmCifBiolAssemblyProvider();
		BioUnitDataProviderFactory.setBioUnitDataProvider(mmcifProvider.getClass().getCanonicalName());	
//		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
	}
	
	private Structure createBioAssembly(PrintWriter error, String pdbId, int i) {
		StructureIO.setAtomCache(cache);
		Structure structure = null;
		try {
			structure = StructureIO.getBiologicalAssembly(pdbId, i);
		} catch (IOException e) {
			// TODO Auto-generated catch block
		//	e.printStackTrace();
			error.println(pdbId + "[" + i + "]: " + e.getMessage());
			error.flush();
		} catch (StructureException e) {
			// TODO Auto-generated catch block
		//	e.printStackTrace();
			error.println(pdbId + "[" + i + "]: " + e.getMessage());
			error.flush();
		}
		return structure;
	}
	
    /** Load a specific biological assembly for a PDB entry
     *  
     * @param pdbId .. the PDB ID
     * @param bioAssemblyId .. the first assembly has the bioAssemblyId 1
     * @return a Structure object or null if something went wrong.
     */
    public Structure readBioAssemblyFile(PrintWriter error, String pdbId, int bioAssemblyId) {
        // pre-computed files use lower case PDB IDs
        pdbId = pdbId.toLowerCase();

        // The low level PDB file parser
        PDBFileReader pdbreader = new PDBFileReader();

        // we just need this to track where to store PDB files
        // this checks the PDB_DIR property (and uses a tmp location if not set) 
        pdbreader.setPath(cache.getPath());

        pdbreader.setFileParsingParameters(cache.getFileParsingParams());

        // download missing files
        pdbreader.setAutoFetch(true);

        pdbreader.setBioAssemblyId(bioAssemblyId);
        pdbreader.setBioAssemblyFallback(false);

        Structure structure = null;
        try { 
            structure = pdbreader.getStructureById(pdbId);
            if ( bioAssemblyId > 0 )
                structure.setBiologicalAssembly(true);
            structure.setPDBCode(pdbId);
        } catch (Exception e){
    //    	e.printStackTrace();
			error.println(pdbId + "[" + bioAssemblyId + "]: " + e.getMessage());
			error.flush();
			return structure;
        }
        return structure;
    }
	
	
//	private static final String[] excludes = new String[]{"1M4X", "2BGJ" , "2J4Z", "2JBP","3HQV","3HR2", "2GSY","2DF7"};
	// 1M4X 540 subunits, WARNING ID 1000> 100, 2081520 atoms, problem with reading BIOMT
	// 2BGJ, small protein, subunits 1, atoms 260 (OK)
	// 2J4Z, small protein, subunits 1, atoms 263 (OK)
	// 2JBP, WARNING ID 3 > 1, 1st BA is single chain: subunits 1, atoms 283
	// 3HQV, 27 subunits, 26865 atoms (OK, local/helical sym. exceeds time limit)
	// 3HR2, 27 subunits, 26757 atoms (sub clusters: 24507, exceeds time limit)
	// 2GSY, 60 subunits, 25680 atoms, made up of 3 AUs, A60/I 7300 ms, AU is asymmetric, exceeds time limit
	// 2DF7, 60 subunits, 24480 atoms, made up of 3 AUs, A60/I 7500 ms
	
	
	
//	private static final String[] excludes = new String[]{"1M4X", "2BGJ" , "2J4Z", "2JBP","3HQV","3HR2", "2GSY","2DF7"};
	private static final String[] testCases = new String[]{"4A1I"};
}
