package org.biojava3.structure.align.symm.quaternary.analysis;

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

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.apache.commons.collections.list.SynchronizedList;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBParseException;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.ChainClusterer;
import org.biojava3.structure.align.symm.quaternary.FindQuarternarySymmetry;
import org.biojava3.structure.align.symm.quaternary.MomentsOfInertia;
import org.biojava3.structure.align.symm.quaternary.QuatSymmetryWriter;
import org.biojava3.structure.align.symm.quaternary.Rotation;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;
import org.biojava3.structure.align.symm.quaternary.Subunits;
import org.biojava3.structure.align.symm.quaternary.MomentsOfInertia.SymmetryClass;
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.rcsb.fatcat.server.PdbChainKey;

import sun.security.acl.GroupImpl;

public class ScanPdbForQuarternarySymmetryNew implements Runnable {
    private int threads = 1;
    private int threadId = 0;
	private static String PDB_PATH = "C:/PDB/";
	private static final int MIN_SEQUENCE_LENGTH = 24;

	public ScanPdbForQuarternarySymmetryNew () {
	}
	
	public ScanPdbForQuarternarySymmetryNew (int threads, int threadId) {
		this.threads = threads;
		this.threadId = threadId;		
		System.out.println("Tread: " + threadId + " threads: " + threads);	
	}

	public void run() {
		AtomCache cache = new AtomCache();
		cache.setAutoFetch(true);
		cache.setPath(PDB_PATH);	
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		p.setLoadChemCompInfo(true);
		//p.setParseCAOnly(true);
		//p.setMaxAtoms(50000000);

		p.setAtomCaThreshold(Integer.MAX_VALUE);
	//	System.out.println("PARSING ALL ATOMS!!!");
		p.setAcceptedAtomNames(new String[]{" CA ", " CB "});
		cache.setFileParsingParams(p);

		System.out.println("Reading blastclust files");

		BlastClustReader reader100 = new BlastClustReader(100);
		BlastClustReader reader90 = new BlastClustReader(90);
		BlastClustReader reader70 = new BlastClustReader(70);
		BlastClustReader reader40 = new BlastClustReader(40);

//		Set<String> reps = GetRepresentatives.getAll();
//		System.out.println("Total representative PDBs: " + reps.size());

		// uncommment the following lines to try just a few examples
//		reps.clear();
		//		reps.add("1AJC");
		// reps.add("1A98"); // case where seq. alignment != structural alignment due to alignment problems at gap.
		//		reps.add("2CD5"); // no bioassembly file available
		//		reps.add("1FT8");
		//		reps.add("4HHB"); // 2 alpha, 2 beta
		//		reps.add("2BG9"); // acetylcholin receptor, 2 alpha, 1 beta, 1 delta, 1 gamma
		//		reps.add("2WRN"); // ribosome
		//		reps.add("3SYW"); // DNA, no protein chains
		//		reps.add("5CSC");
		//		reps.add("1B77");
		//		reps.add("1R0B");
		//		reps.add("2FZC");
		//		reps.add("4HHB");

		//		reps.add("1JRE"); // (1JRE_A)12: A12:T
		//		reps.add("1JTS");
		//		reps.add("1L8H");
		//		reps.add("1L8I");

		//		reps.add("4HHB");
		//		reps.add("3R5I");

		//		reps.add("1AA1");//  A8B8:D4
		//		reps.add("1AUS");
		//		reps.add("1RBO");
		//		reps.add("1RCO");
		//		reps.add("1RCX");
		//		reps.add("1RXO");
		//		reps.add("1UPM");
		//		reps.add("1UPP");
		//		reps.add("8RUC");
		//		
		//		reps.add("1AHU"); // A8:D4
		//		reps.add("1IXO");
		//		reps.add("2BFO");
		//		reps.add("3H4V");
		//		reps.add("2BFM");
		//		reps.add("1W0C");

//	    boolean writeFile = true;
		boolean writeFile = false;
		PrintWriter out = null;
		PrintWriter out1 = null;
		PrintWriter error = null;
		try {
			out = new PrintWriter(new FileWriter(PDB_PATH + "rep_sym" + threadId + ".csv"));
			out1 = new PrintWriter(new FileWriter(PDB_PATH + "rep_err" + threadId + ".csv"));
			error = new PrintWriter(new FileWriter(PDB_PATH + "error" + threadId + ".txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}


		long t1 = System.nanoTime();
		long symTime = 0;

		int multimer = 0;
		int excluded = 0;
		int total = 0;
		int err = 0;

		out.println("pdbId,formula,signature100,stoichiometry100,types100,signature90,stoichiometry90,types90,signature70,stoichiometry70,types70,signature40,stoichiometry40,types40,pointgroup,symops,multiplicity,cacount,chains,method,rmsdS,rmsdT,gts,symclass,asymcoeff,bioassembly,time,jmol,trace");
		out1.println("pdbId,formula,signature100,stoichiometry100,types100,signature90,stoichiometry90,types90,signature70,stoichiometry70,types70,signature40,stoichiometry40,types40,pointgroup,complete,symops,multiplicity,cacount,chains,method,rmsdS,rmsdT,gts,symclass,asymcoeff,bioassembly,time,seqNum,unksequence,jmol,trace");

		System.out.println("Getting PdbEntryInfo");
		List<PdbEntryInfo> list = PdbEntryInfoParser.getPdbEntryInfo();
        int chunkSize = (int) Math.round(Math.ceil((double)list.size()/threads));
        int start = threadId * chunkSize;
        int end = start + chunkSize - 1;
        end = Math.min(end, list.size());
  //      end = start + 20;
        System.out.println("Start: " + start);
        System.out.println("End: " + end);
        
        boolean skip = false;
        String restartId = "3THS";

		for (int k = start; k < end; k++) {	
			PdbEntryInfo entry = list.get(k);
			total++;
			String pdbId = entry.getPdbId();
			
			if (skip && pdbId.equals(restartId)) {
				skip = false;
			} 
			
			if (skip) {
				continue;
			}
//			if (!pdbId.equals("4HHB")) continue; 
			if (!pdbId.equals("2BG9")) continue; 
			// 1B4N, 1AVO, 1BZ5, 1A8S, 1BHC(D5), 1KN0(C11) // good example
//			if (!pdbId.equals("1OCW")) continue; 
//			if (!pdbId.equals("1A6D")) continue; // trace <= 0
//			if (!pdbId.equals("1A17")) continue;
//			if (!pdbId.equals("4HHB")) continue; // 
//			if (!pdbId.equals("3KO1")) continue; // good
//			if (!pdbId.equals("2AO9")) continue; // good
//			if (!pdbId.equals("3C9K")) continue; 
//			if (!pdbId.equals("1A30")) continue;
//			if (!pdbId.equals("2CB2")) continue;
//			if (!pdbId.equals("1AEW")) continue;
//			if (!pdbId.equals("1NF4")) continue; // ???? octahedral//
//			if (!pdbId.equals("1H2I")) continue; // good
//			if (!pdbId.equals("1A34")) continue;
//			if (!pdbId.equals("1P3H")) continue;
//			if (!pdbId.equals("1M5Q")) continue; // good example
//			if (!pdbId.equals("1A4Q")) continue;
//			if (!pdbId.equals("1M4X")) continue;
//			if (!pdbId.equals("1ML5")) continue;
//			if (!pdbId.equals("3LJ4")) continue;
			System.out.println("------------- " + pdbId  + "-------------");

			int bioAssemblyCount = entry.getBioAssemblyCount();
			System.out.println("Bioassemblies: " + bioAssemblyCount);
			int n = Math.max(bioAssemblyCount, 1);
			for (int i = 0; i < n; i++) {
	//			if (pdbId.equals("1M4X")) continue; // largest PDB assembly, causes occasional GC error

				Structure structure = null;
				try {
					if (bioAssemblyCount == 0) {
						structure = cache.getStructure(pdbId);
					} else {
						structure = cache.getBiologicalAssembly(pdbId, i+1, true);
					}

				} catch (StructureException e) {
					e.printStackTrace();
					error.println(pdbId + "------------------------------------");
					error.println(e.getMessage());
					error.flush();
					err++;
					continue;
				} catch (IOException e) {
					e.printStackTrace();
					error.println(pdbId + "------------------------------------");
					error.println(e.getMessage());
					error.flush();
					err++;
					continue;
				}

				boolean biologicalAssembly = structure.isBiologicalAssembly();

				long tc1 = System.nanoTime(); 	

				System.out.println("Loaded 2bg9");
				FindQuarternarySymmetry finder = new FindQuarternarySymmetry(structure);
				finder.setMinimumSequenceLength(MIN_SEQUENCE_LENGTH);
				
				finder.setPseudoSymmetryAllowed(false);
				RotationGroup rotationGroup = finder.getRotationGroup();	
				String pointGroup = rotationGroup.getPointGroup();
				System.out.println("Point group: " + pointGroup);
				
				ChainClusterer grouper = finder.getSequenceCluster();
				String formula = grouper.getCompositionFormula();
				System.out.println("Formula: " + formula);
				boolean unknownSequence = grouper.containsUnknownSequence();
				System.out.println(grouper);

				String method = finder.getMethod();	

				// get metrics
				int caCount = finder.getSubunits().getCalphaCount();
				boolean groupComplete = rotationGroup.isComplete();
				float rmsd = (float) rotationGroup.getAverageSubunitRmsd();
				float rmsdT = (float) rotationGroup.getAverageTraceRmsd();
				float gts = (float) rotationGroup.getAverageTraceGtsMin();
		
				int order = rotationGroup.getOrder();
				float multiplicity = 1;

				if (rotationGroup.getOrder() > 0) {
					multiplicity = grouper.getMultiplicity()/(float)order;
				}
				
				Subunits subunits = finder.getSubunits();
				Matrix4d matrix = AxisTransformation.getTransformation(structure, subunits, rotationGroup);
				System.out.println("Transformation:");
	//			System.out.println(matrix);
				String jmol = AxisTransformation.getJmolQuat(matrix);
				System.out.println(jmol);
				double trace = AxisTransformation.getTrace(matrix);
				
				// determine overall symmetry
			
				MomentsOfInertia m = subunits.getMomentsOfInertia();
				MomentsOfInertia.SymmetryClass symmetryClass = m.getSymmetryClass(0.2);
				double asymmetryCoefficient = m.getAsymmetryParameter(0.001);

				// create list of representative chains
				ProteinComplexSignature s100 = new ProteinComplexSignature(pdbId, grouper, reader100);
				String signature100 = s100.getComplexSignature();
				String stoich100 = s100.getComplexStoichiometry();
				int types100 = s100.getSubunitTypeCount();
				ProteinComplexSignature s90 = new ProteinComplexSignature(pdbId, grouper, reader90);
				String signature90 = s90.getComplexSignature();
				String stoich90 = s90.getComplexStoichiometry();
				int types90 = s90.getSubunitTypeCount();
				System.out.println("sign100: " + signature100 + ":" + stoich100);
				ProteinComplexSignature s70 = new ProteinComplexSignature(pdbId, grouper, reader70);
				String signature70 = s70.getComplexSignature();
				String stoich70 = s70.getComplexStoichiometry();
				int types70 = s70.getSubunitTypeCount();
				System.out.println("sign70: " + signature70 + ":" + stoich70);
				ProteinComplexSignature s40 = new ProteinComplexSignature(pdbId, grouper, reader40);
				String signature40 = s40.getComplexSignature();
				String stoich40 = s40.getComplexStoichiometry();
				int types40 = s40.getSubunitTypeCount();
				System.out.println("sign40: " + signature40 + ":" + stoich40);
				
				// TODO use chain signatures to label interactions of ligands
				long tc2 = System.nanoTime();
				long time = (tc2 - tc1)/1000000;
				symTime += time;


				
				// write .csv summary file
				if (finder.getChainCount() > 1 && groupComplete && !unknownSequence) {			
					out.print(pdbId + "," + formula + "," + signature100 + "," + stoich100 + "," + types100 + "," + signature90 + "," + stoich90 + "," + types90 + "," + signature70 + "," + stoich70  + "," + types70 + "," + signature40 + "," + stoich40 + "," + types40 + "," + pointGroup + "," +
							order + "," + multiplicity + "," + caCount + "," + finder.getChainCount() + "," + method  + "," + rmsd + "," + rmsdT + "," + gts + "," + symmetryClass + "," + asymmetryCoefficient + "," +
							biologicalAssembly + "," + time + "," + jmol + "," + trace);
							out.println();
							out.flush();
							multimer++;
				} else if (finder.getChainCount() > 1) {
					out1.print(pdbId + "," + formula + "," + signature100 + "," + stoich100 + "," + types100 + "," + signature90 + "," + stoich90 + "," + types90 + "," + signature70 + "," + stoich70  + "," + types70 + "," + signature40 + "," + stoich40 + "," + types40 + "," + pointGroup + "," + groupComplete + "," +
							order + "," + multiplicity + "," + caCount + "," + finder.getChainCount() + "," + "," + method  + "," + rmsd + "," + rmsdT + "," + gts + "," + symmetryClass + "," + asymmetryCoefficient + "," +
							biologicalAssembly + "," + time + "," + "," + unknownSequence + "," + jmol + "," + trace);
					out1.println();
					out1.flush();
					excluded++;
				}

				//			if (total == 500) break;		
				
				// write symmetry copies to file
				if (writeFile) {
					QuatSymmetryWriter writer = new QuatSymmetryWriter(structure);
					writer.writeTransformedStructure(matrix, "C:/PDB/" + pdbId + "_trans.pdb");
					structure = new StructureImpl();
					addAxis(structure, subunits, rotationGroup);
					QuatSymmetryWriter writer1 = new QuatSymmetryWriter(structure);
	//				PermutationGroup pgroup = new PermutationGroup();
					System.out.println("Rotation group order: " + rotationGroup.getOrder());
					for (int j = 0; j < rotationGroup.getOrder();j++) {
						//					System.out.println("Matrix:");
						String fileName = PDB_PATH + "transformed/" + pdbId + j + ".pdb";
						//					System.out.println("Direction:");
						//					System.out.println(rotationGroup.getRotation(j).getDirection());
						//					System.out.println("Fold:");
						//					System.out.println(rotationGroup.getRotation(j).getFold());
						//					System.out.println("Permutation:");
						//					System.out.println(rotationGroup.getRotation(j).getPermutation());
						//					pgroup.addPermutation(rotationGroup.getRotation(j).getPermutation());
						//			System.out.println("Transformation:");
						//			System.out.println(rotationGroup.getRotation(j).getTransformation());
						//					System.out.println("Axis:");
						//					System.out.println(Math.toDegrees(rotationGroup.getRotation(j).getAxisAngle().angle));
						//					System.out.println(rotationGroup.getRotation(j).getAxisAngle());
						//					System.out.println("writing transformed file: " + fileName);
											writer1.writeTransformedStructure(rotationGroup.getRotation(j).getTransformation(), fileName);	
					}
				}
				
			}
		}
		long t2 = System.nanoTime();

		System.out.println("Cpu time: " + (t2-t1)/1000000 + " ms.");
		System.out.println("Calc. time: " + symTime);
		System.out.println("Multimers: " + multimer + "out of " + total);
		System.out.println("Excluded: " + excluded);
		System.out.println("Errors: " + err);
		System.out.println("Total structure: " + PdbEntryInfoParser.getPdbEntryInfo().size());
		out.close();
		out1.close();
		error.close();
	}

	private static void addAxis(Structure structure, Subunits subunits, RotationGroup rotGroup) {
		Group g = new HetatomImpl();
		try {
			g.setPDBName("AXS");
		} catch (PDBParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		g.setResidueNumber("X", 1, ' ');

		Chain c = new ChainImpl();
		c.addGroup(g);
		c.setChainID("X");	
		structure.addChain(c);

		int count = 1;

		for (int i = 0; i < rotGroup.getOrder(); i++) {
			Rotation rotation = rotGroup.getRotation(i);
			AxisAngle4d axisAngle1 = rotation.getAxisAngle();
			Vector3d axisVector = new Vector3d(axisAngle1.x, axisAngle1.y, axisAngle1.z);
			Point3d center = subunits.getCentroid();

			Vector3d axisVector1 = new Vector3d(axisVector);
			axisVector1.scaleAdd(25.0, center);
			Atom dummy1 = new AtomImpl();
			double[] c1 = {axisVector1.x, axisVector1.y, axisVector1.z};
			dummy1.setCoords(c1);
			dummy1.setElement(Element.Xe);
			dummy1.setFullName("XE1 ");
			dummy1.setPDBserial(count);
			count++;
			g.addAtom(dummy1);

			Vector3d axisVector2 = new Vector3d(axisVector);
			axisVector2.scaleAdd(-25.0, center);
			Atom dummy2 = new AtomImpl();
			double[] c2 = {axisVector2.x, axisVector2.y, axisVector2.z};
			dummy2.setCoords(c2);
			dummy2.setElement(Element.Xe);
			dummy2.setFullName("XE2 ");
			dummy2.setPDBserial(count);
			count++;
			g.addAtom(dummy2);			
		}
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
