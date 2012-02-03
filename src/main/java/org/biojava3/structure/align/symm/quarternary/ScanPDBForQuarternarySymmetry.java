package org.biojava3.structure.align.symm.quarternary;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
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

public class ScanPDBForQuarternarySymmetry {

	public static String PDB_PATH = "C:/PDB/";
	
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
		p.setAcceptedAtomNames(new String[]{" CA ", " CB "});
		cache.setFileParsingParams(p);

		Map<Integer,String> pointGroupMap = new HashMap<Integer,String>();


		Set<String> reps = GetRepresentatives.getAll();
		System.out.println("Total representative PDBs: " + reps.size());

		// uncommment the following lines to try just a few examples
//		reps.clear();
//		reps.add("1WG8");

		boolean writeFile = false;
		PrintWriter out = null;
		PrintWriter error = null;
		try {
			out = new PrintWriter(new FileWriter(PDB_PATH + "rep_sym.csv"));
			error = new PrintWriter(new FileWriter(PDB_PATH + "error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}


		long t1 = System.nanoTime();
		GlobalSequenceGrouper grouper = new GlobalSequenceGrouper();

		int multimer = 0;
		int total = 0;
		for (String pdbId : reps){
			System.out.println("------------- " + pdbId  + "-------------");
			if (total == 0) {
				out.println("pdbId,formula,cacount,chains,hashcode,match,method,pointgroup,symops,multiplicity,rmsdS,rmsdT,gts,symclass,asymcoeff,seqclusters,pgcomplete,bioassembly,time,seqNum");
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
			grouper.setStructure(structure);
			int minSequenceLength = 24;
			grouper.setMinSequenceLength(minSequenceLength);
			int nClusters = grouper.getSequenceCluster100().size();
			String formula = grouper.getCompositionFormula();
			int hashCode = grouper.hashCodeMD5();
			boolean sequenceNumberedCorrectly = grouper.isSequenceNumberedCorrectly();
	
			// determine point group
			List<Point3d[]> caCoords = grouper.getCalphaTraces();
			List<Point3d[]> cbCoords = grouper.getCbetaTraces();
			List<Integer> sequenceClusterIds = grouper.getSequenceClusterIds();
			FindQuarternarySymmetry finder = new FindQuarternarySymmetry(caCoords, cbCoords, sequenceClusterIds);
			finder.setPseudoSymmetryAllowed(false);
			RotationGroup rotationGroup = finder.getRotationGroup();	
			String pointGroup = rotationGroup.getPointGroup();
			
			// check if all complexes with the same composition have the same point group
			String pg = pointGroupMap.get(hashCode);
			boolean consistent = true;
			if (pg == null) {
				pointGroupMap.put(hashCode, pointGroup);
			} else {
				if (!pg.equals(pointGroup)) {
					consistent = false;
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
			double multiplicity = 1;
			if (rotationGroup.getOrder() > 0) {
				multiplicity = finder.getChainCount()/rotationGroup.getOrder();
			}

			long tc2 = System.nanoTime();
			long time = (tc2 - tc1)/1000000;
			
			// write .csv summary file
			if (finder.getChainCount() > 1) {			
				out.println(pdbId + "," + formula + "," + caCount + "," + finder.getChainCount() + "," + hashCode +  "," + consistent + "," + method + "," + pointGroup + "," + 
						rotationGroup.getOrder()+ "," + multiplicity + "," + rmsd + "," + rmsdT + "," + gts + "," + symmetryClass + "," + asymmetryCoefficient + "," +
						nClusters + "," + groupComplete + "," + 
						structure.isBiologicalAssembly() + "," + time + "," + sequenceNumberedCorrectly);
				out.flush();
				multimer++;
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
