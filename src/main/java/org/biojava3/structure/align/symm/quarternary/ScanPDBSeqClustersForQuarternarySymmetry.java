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

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

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
import org.biojava3.structure.dbscan.GetRepresentatives;
import org.rcsb.fatcat.server.PdbChainKey;

import sun.security.acl.GroupImpl;

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

        BlastClustReader reader100 = new BlastClustReader(100);
        BlastClustReader reader90 = new BlastClustReader(90);
        BlastClustReader reader70 = new BlastClustReader(70);
        BlastClustReader reader40 = new BlastClustReader(40);
        
		Set<String> reps = GetRepresentatives.getAll();
		System.out.println("Total representative PDBs: " + reps.size());

		// uncommment the following lines to try just a few examples
//		reps.clear();
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

//		boolean writeFile = false;
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
		long symTime = 0;

		int multimer = 0;
		int excluded = 0;
		int total = 0;
		int err = 0;
		
		
		for (String pdbId : reps){
			System.out.println("------------- " + pdbId  + "-------------");
			if (total == 0) {
				// TODO number of nucleic acid contacts?
				out.println("pdbId,formula,signature100,stoichiometry100,types100,signature90,stoichiometry90,types90,signature70,stoichiometry70,types70,signature40,stoichiometry40,types40,pointgroup,symops,multiplicity,isoquaternary,cacount,chains,method,rmsdS,rmsdT,gts,symclass,asymcoeff,bioassembly,time,ligands,interactiontype");
				out1.println("pdbId,formula,signature100,stoichiometry100,types100,signature90,stoichiometry90,types90,signature70,stoichiometry70,types70,signature40,stoichiometry40,types40,pointgroup,complete,symops,multiplicity,isoquaternary,cacount,chains,method,rmsdS,rmsdT,gts,symclass,asymcoeff,bioassembly,time,seqNum,unksequence,ligands,interactiontype");
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

			// cluster sequences by sequence identity
			GlobalSequenceGrouper grouper = new GlobalSequenceGrouper(structure, MIN_SEQUENCE_LENGTH);
			String formula = grouper.getCompositionFormula();
			System.out.println("Formula: " + formula);
			boolean sequenceNumberedCorrectly = grouper.isSequenceNumberedCorrectly();
			boolean unknownSequence = grouper.isUnknownSequence();
			
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
	
			// determine point group
			List<Point3d[]> caCoords = grouper.getCalphaCoordinates();
			List<Point3d[]> cbCoords = grouper.getCbetaCoordinates();
			List<Integer> sequenceClusterIds = grouper.getSequenceClusterIds();
			FindQuarternarySymmetry finder = new FindQuarternarySymmetry(caCoords, cbCoords, sequenceClusterIds);
			finder.setPseudoSymmetryAllowed(false);
			RotationGroup rotationGroup = finder.getRotationGroup();	
			String pointGroup = rotationGroup.getPointGroup();
			
			LigandInteractions li = new LigandInteractions(structure, s100);
			li.setInteractingChains(grouper.getChains());
			
			// check if all complexes with the same composition have the same point group
			String pg = pointGroupMap.get(signature100);
			boolean isoQuaternary = true;
			if (pg == null) {
				pointGroupMap.put(signature100, pointGroup);
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
			MomentsOfInertia.SymmetryClass symmetryClass = m.getSymmetryClass(0.2);
			double asymmetryCoefficient = m.getAsymmetryParameter(0.001);
			
			int order = rotationGroup.getOrder();
			float multiplicity = 1;
	
			if (rotationGroup.getOrder() > 0) {
				multiplicity = grouper.getMultiplicity()/(float)order;
			}

			long tc2 = System.nanoTime();
			long time = (tc2 - tc1)/1000000;
			symTime += time;
			
			// write .csv summary file
			if (finder.getChainCount() > 1 && groupComplete && sequenceNumberedCorrectly && !unknownSequence) {			
				out.print(pdbId + "," + formula + "," + signature100 + "," + stoich100 + "," + types100 + "," + signature90 + "," + stoich90 + "," + types90 + "," + signature70 + "," + stoich70  + "," + types70 + "," + signature40 + "," + stoich40 + "," + types40 + "," + pointGroup + "," +
						order + "," + multiplicity + "," + isoQuaternary + "," + caCount + "," + finder.getChainCount() + "," + method  + "," + rmsd + "," + rmsdT + "," + gts + "," + symmetryClass + "," + asymmetryCoefficient + "," +
						biologicalAssembly + "," + time);
				out.print("," + li.toString());
				out.print("," + li.getInteractionType());
				out.println();
				out.flush();
				multimer++;
			} else if (finder.getChainCount() > 1) {
				out1.print(pdbId + "," + formula + "," + signature100 + "," + stoich100 + "," + types100 + "," + signature90 + "," + stoich90 + "," + types90 + "," + signature70 + "," + stoich70  + "," + types70 + "," + signature40 + "," + stoich40 + "," + types40 + "," + pointGroup + "," + groupComplete + "," +
						order + "," + multiplicity + "," + isoQuaternary + "," + caCount + "," + finder.getChainCount() + "," + "," + method  + "," + rmsd + "," + rmsdT + "," + gts + "," + symmetryClass + "," + asymmetryCoefficient + "," +
						biologicalAssembly + "," + time + "," + sequenceNumberedCorrectly + "," + unknownSequence);
				out1.print("," + li.toString());
				out1.print("," + li.getInteractionType());
				out1.println();
				out1.flush();
				excluded++;
			}

			//			if (total == 500) break;

			// write symmetry copies to file
			if (writeFile) {
				structure = new StructureImpl();
				addAxis(structure, subunits, rotationGroup);
				QuatSymmetryWriter writer = new QuatSymmetryWriter(structure);
				PermutationGroup pgroup = new PermutationGroup();
				System.out.println("Rotation group order: " + rotationGroup.getOrder());
				for (int i = 0; i < rotationGroup.getOrder(); i++) {
//					System.out.println("Matrix:");
					String fileName = PDB_PATH + "transformed/" + pdbId + i + ".pdb";
//					System.out.println("Direction:");
//					System.out.println(rotationGroup.getRotation(i).getDirection());
//					System.out.println("Fold:");
//					System.out.println(rotationGroup.getRotation(i).getFold());
//					System.out.println("Permutation:");
//					System.out.println(rotationGroup.getRotation(i).getPermutation());
//					pgroup.addPermutation(rotationGroup.getRotation(i).getPermutation());
		//			System.out.println("Transformation:");
		//			System.out.println(rotationGroup.getRotation(i).getTransformation());
//					System.out.println("Axis:");
//					System.out.println(Math.toDegrees(rotationGroup.getRotation(i).getAxisAngle().angle));
//					System.out.println(rotationGroup.getRotation(i).getAxisAngle());
//					System.out.println("writing transformed file: " + fileName);
//					writer.writeTransformedStructure(rotationGroup.getRotation(i).getTransformation(), fileName);	
				}
//				System.out.println("Group hashcode: " + pgroup.hashCode());
//				System.out.println(pgroup.getGroupTable());
			}
			total++;
		}
		long t2 = System.nanoTime();
		
		System.out.println("Cpu time: " + (t2-t1)/1000000 + " ms.");
		System.out.println("Calc. time: " + symTime);
		System.out.println("Multimers: " + multimer + "out of " + total);
		System.out.println("Excluded: " + excluded);
		System.out.println("Errors: " + err);
		System.out.println("Total structure: " + reps.size());
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
