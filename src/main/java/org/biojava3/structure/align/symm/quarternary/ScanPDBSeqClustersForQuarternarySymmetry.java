package org.biojava3.structure.align.symm.quarternary;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.rcsb.fatcat.server.PdbChainKey;

public class ScanPDBSeqClustersForQuarternarySymmetry {

	private static String PDB_PATH = "C:/PDB/";
	private static final int MIN_SEQUENCE_LENGTH = 24;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		AtomCache cache = new AtomCache();
		cache.setAutoFetch(true);
		cache.setPath(PDB_PATH);	
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		//p.setParseCAOnly(true);
		//p.setMaxAtoms(50000000);
		
		p.setAtomCaThreshold(Integer.MAX_VALUE);
		System.out.println("PARSING ALL ATOMS!!!");
	//	p.setAcceptedAtomNames(new String[]{" CA ", " CB "});
		cache.setFileParsingParams(p);

		Map<String,String> pointGroupMap = new HashMap<String,String>();

        BlastClustReader reader = new BlastClustReader(100);
        
		Set<String> reps = GetRepresentatives.getAll();
		System.out.println("Total representative PDBs: " + reps.size());

		// uncommment the following lines to try just a few examples
//		reps.clear();
//		reps.add("1AA7");
//		reps.add("1B4A");
//		reps.add("1C4U");
//		reps.add("3IYN");
//		reps.add("4HHB");
//		reps.add("1IRD");
//		reps.add("2DN2");
//		reps.add("1E7W");
//		reps.add("1E92");
//		reps.add("2BF7");
//		reps.add("2BFP");
//		reps.add("2BFA");
//		reps.add("2BFO");
//		reps.add("3H4V");
//		reps.add("2BFM");
//		reps.add("1W0C");

		boolean writeFile = false;
		PrintWriter out = null;
		PrintWriter out1 = null;
		PrintWriter error = null;
		try {
			out = new PrintWriter(new FileWriter(PDB_PATH + "rep_sym.csv"));
			out1 = new PrintWriter(new FileWriter(PDB_PATH + "rep_err.csv"));
			error = new PrintWriter(new FileWriter(PDB_PATH + "error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}


		long t1 = System.nanoTime();

		int multimer = 0;
		int excluded = 0;
		int total = 0;
		
		
		for (String pdbId : reps){
			System.out.println("------------- " + pdbId  + "-------------");
			if (total == 0) {
				out.println("pdbId,formula,signature,pointgroup,symops,multiplicity,isoquaternary,cacount,chains,method,rmsdS,rmsdT,gts,symclass,asymcoeff,bioassembly,time,ligands");
				out1.println("pdbId,formula,signature,pointgroup,complete,symops,multiplicity,isoquaternary,cacount,chains,method,rmsdS,rmsdT,gts,symclass,asymcoeff,bioassembly,time,seqNum,unksequence,ligands");
			}

			if (pdbId.equals("1M4X")) continue; // largest PDB assembly, causes occasional GC error

			Structure structure = null;
			try {
				structure = cache.getBiologicalAssembly(pdbId, 1, true);
				//				structure = cache. getStructure(pdbId);

			} catch (StructureException e) {
				e.printStackTrace();
				error.println(pdbId + "------------------------------------");
				error.println(e.getMessage());
				error.flush();
				continue;
			} catch (IOException e) {
				e.printStackTrace();
				error.println(pdbId + "------------------------------------");
				error.println(e.getMessage());
				error.flush();
				continue;
			}
			
			boolean biologicalAssembly = structure.isBiologicalAssembly();

			long tc1 = System.nanoTime();

			// cluster sequences by sequence identity
			GlobalSequenceGrouper grouper = new GlobalSequenceGrouper(structure, MIN_SEQUENCE_LENGTH);
			String formula = grouper.getCompositionFormula();
			System.out.println("Formula: " + formula);
			boolean sequenceNumberedCorrectly = grouper.isSequenceNumberedCorrectly();
			boolean unknownSequence = grouper.isUnknownSequence();
			
			// create list of representative chains
			ProteinComplexSignature s = new ProteinComplexSignature(pdbId, grouper, reader);
            String signature = s.getComplexSignature();
            // TODO use chain signatures to label interactions of ligands
	
			// determine point group
			List<Point3d[]> caCoords = grouper.getCalphaCoordinates();
			List<Point3d[]> cbCoords = grouper.getCbetaCoordinates();
			List<Integer> sequenceClusterIds = grouper.getSequenceClusterIds();
			FindQuarternarySymmetry finder = new FindQuarternarySymmetry(caCoords, cbCoords, sequenceClusterIds);
			finder.setPseudoSymmetryAllowed(false);
			RotationGroup rotationGroup = finder.getRotationGroup();	
			String pointGroup = rotationGroup.getPointGroup();
			
			LigandInteractions li = new LigandInteractions(structure, s);
			li.setInteractingChains(grouper.getChains());
//			List<InteractingLigand> ligands = li.getInteractingLigands();
//			List<String> iLigs = new ArrayList<String>();
//			for (InteractingLigand lig: ligands) {
//				iLigs.add(lig.toString());
//			}
//			Collections.sort(iLigs);
			
			
			// check if all complexes with the same composition have the same point group
			String pg = pointGroupMap.get(signature);
			boolean isoQuaternary = true;
			if (pg == null) {
				pointGroupMap.put(signature, pointGroup);
			} else {
				if (!pg.equals(pointGroup)) {
					isoQuaternary = false;
				}
			}
			String method = finder.getMethod();	
			
			// get metrics
			int caCount = finder.getSubunits().getCalphaCount();
			boolean groupComplete = rotationGroup.isComplete();
			float rmsd = (float) rotationGroup.getAverageSubunitRmsd();
			float rmsdT = (float) rotationGroup.getAverageTraceRmsd();
			float gts = (float) rotationGroup.getAverageTraceGtsMin();
			
			
			// determine overall symmetry
			Subunits subunits = finder.getSubunits();
			MomentsOfInertia m = subunits.getMomentsOfInertia();
			MomentsOfInertia.SymmetryClass symmetryClass = m.getSymmetryClass(0.05);
			double asymmetryCoefficient = m.getAsymmetryParameter(0.05);
			
			int order = rotationGroup.getOrder();
			float multiplicity = 1;
	
			if (rotationGroup.getOrder() > 0) {
				multiplicity = grouper.getMultiplicity()/(float)order;
			}

			long tc2 = System.nanoTime();
			long time = (tc2 - tc1)/1000000;
			
			// write .csv summary file
			if (finder.getChainCount() > 1 && groupComplete && sequenceNumberedCorrectly && !unknownSequence) {			
				out.print(pdbId + "," + formula + "," + signature + "," + pointGroup + "," +
						order + "," + multiplicity + "," + isoQuaternary + "," + caCount + "," + finder.getChainCount()  + "," + method  + "," + rmsd + "," + rmsdT + "," + gts + "," + symmetryClass + "," + asymmetryCoefficient + "," +
						biologicalAssembly + "," + time);
		//		for (String sl: iLigs) {
		//			out.print("," + sl);
		//		}
				out.print("," + li.toString());
				out.println();
				out.flush();
				multimer++;
			} else if (finder.getChainCount() > 1) {
				out1.print(pdbId + "," + formula + "," + signature + "," + pointGroup + "," + groupComplete + "," +
						order + "," + multiplicity + "," + isoQuaternary + "," + caCount + "," + finder.getChainCount()  + "," + method  + "," + rmsd + "," + rmsdT + "," + gts + "," + symmetryClass + "," + asymmetryCoefficient + "," +
						biologicalAssembly + "," + time + "," + sequenceNumberedCorrectly + "," + unknownSequence);
		//		for (String sl: iLigs) {
		//			out.print("," + sl);
		//		}
				out1.print("," + li.toString());
				out1.println();
				out1.flush();
				excluded++;
			}

			//			if (total == 500) break;

			// write symmetry copies to file
			if (writeFile) {
				QuatSymmetryWriter writer = new QuatSymmetryWriter(structure);
				for (int i = 0; i < rotationGroup.getOrder(); i++) {
					String fileName = PDB_PATH + "transformed/" + pdbId + i + ".pdb";
					writer.writeTransformedStructure(rotationGroup.getRotation(i).getTransformation(), fileName);	
				}
			}
			total++;
		}
		long t2 = System.nanoTime();
		
		System.out.println("Cpu Time: " + (t2-t1)/1000000 + " ms.");
		System.out.println("Multimers: " + multimer + "out of " + total);
		System.out.println("Total structure: " + reps.size());
		out.close();
		out1.close();
		error.close();
	}

	private static Set<String> getRepresentativePDBIds(int sequenceIdentity) {
		System.out.println("GetRepresentatives");
		Set<String> set = new TreeSet<String>();
		SortedSet<PdbChainKey> reps = GetRepresentatives.getRepresentatives(sequenceIdentity);
		for ( PdbChainKey r : reps){
			set.add(r.getPdbId());
		}
		return set;
	}

	private static Set<String> getProblemCases() {
		Set<String> set = new LinkedHashSet<String>();
		set.addAll(Arrays.asList("2W49", "3LUE", "3IKU", "3HQV", "3EDL", "3BYH", "2OM3", "2V6L", "2HIL", 
				"2C0W", "1HGV", "1QL1", "4IFM", "1CGM", "2WYY"));
		return set;
	}

	private static Set<String> getHelicalStructures() {
		Set<String> set = new LinkedHashSet<String>();
		set.addAll(Arrays.asList("2W49", "3LUE", "3IKU", "3HQV", "3EDL", "3BYH", "2OM3", "2V6L", "2HIL", 
				"2C0W", "1HGV", "1QL1", "4IFM", "1CGM", "2WYY"));
		return set;
	}

}
