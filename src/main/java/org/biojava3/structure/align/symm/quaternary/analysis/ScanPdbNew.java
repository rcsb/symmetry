package org.biojava3.structure.align.symm.quaternary.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Set;

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
import org.biojava3.structure.align.symm.geometry.MomentsOfInertia;
import org.biojava3.structure.align.symm.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava3.structure.align.symm.quaternary.AxisTransformation;
import org.biojava3.structure.align.symm.quaternary.FindQuarternarySymmetry;
import org.biojava3.structure.align.symm.quaternary.QuatSymmetryParameters;
import org.biojava3.structure.align.symm.quaternary.QuatSymmetryWriter;
import org.biojava3.structure.align.symm.quaternary.Rotation;
import org.biojava3.structure.align.symm.quaternary.RotationGroup;
import org.biojava3.structure.align.symm.quaternary.Subunits;
import org.biojava3.structure.dbscan.GetRepresentatives;

public class ScanPdbNew implements Runnable {
	private static String PDB_PATH = "C:/Users/Peter/Documents/PDB/";
	private static final int MIN_CHAINS = 1;
	private static final int MIN_SEQUENCE_LENGTH = 24;
	private static final double SEQUENCE_IDENTITY_THRESHOLD = 0.30;
	private static final double ALIGNMENT_FRACTION_THRESHOLD = 0.9;
	private static final double RMSD_THRESHOLD = 5.0;

	public ScanPdbNew () {
	}
	
	public static void main(String[] args) {
		new ScanPdbNew().run();
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
		BlastClustReader reader95 = new BlastClustReader(95);
		BlastClustReader reader70 = new BlastClustReader(70);
		BlastClustReader reader40 = new BlastClustReader(40);

		//		reps.add("1AJC");
		// reps.add("1A98"); // case where seq. alignment != structural alignment due to alignment problems at gap.
		//		reps.add("2CD5"); // no bioassembly file available
		//		reps.add("4HHB"); // 2 alpha, 2 beta
		//		reps.add("2BG9"); // acetylcholin receptor, 2 alpha, 1 beta, 1 delta, 1 gamma
		//		reps.add("2WRN"); // ribosome
		//		reps.add("3SYW"); // DNA, no protein chains
		
	

//	    boolean writeFile = true;
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

		out.println("pdbId,bioassembly,formula,signature100,stoichiometry100,types100,signature95,stoichiometry95,types95,signature70,stoichiometry70,types70,signature40,stoichiometry40,types40,pointgroup,symops,cacount,chains,rmsdS,rmsdT,symclass,asymcoeff,time,rx,ry,rz,bin2,jmol,jmolaxes");
		out1.println("pdbId,bioassembly,formula,signature100,stoichiometry100,types100,signature95,stoichiometry95,types95,signature70,stoichiometry70,types70,signature40,stoichiometry40,types40,pointgroup,complete,symops,cacount,chains,rmsdS,rmsdT,symclass,asymcoeff,time,rx,ry,rz,bin2,jmol,jmolaxes");

		System.out.println("Getting PdbEntryInfo");
//		List<PdbEntryInfo> list = PdbEntryInfoParser.getPdbEntryInfo();
      
        boolean skip = false;
        String restartId = "3THS";
        Set<String> set = GetRepresentatives.getAll();

//		for (int k = 0; k < list.size(); k++) {	
//			PdbEntryInfo entry = list.get(k);
//			total++;
//			String pdbId = entry.getPdbId();
        for (String pdbId: set) {
			
			if (skip && pdbId.equals(restartId)) {
				skip = false;
			} 
			
			if (skip) {
				continue;
			}

//			if (!pdbId.equals("3LSV")) continue; // C1 (A)	
//			if (!pdbId.equals("1VFB")) continue; // C1 (A)	
//			if (!pdbId.equals("1J1X")) continue; // C1 (A)
//			if (!pdbId.equals("1B27")) continue; // C2(AB) 2 BAs
//			if (!pdbId.equals("1S6V")) continue; // C1(AB) 2 BAs
//			if (!pdbId.equals("1AFA")) continue;
//			if (!pdbId.equals("3T88")) continue;
//			if (!pdbId.equals("1A9S")) continue; // good example
//			if (!pdbId.equals("1B44")) continue;
//			if (!pdbId.equals("4EAM")) continue;
//			if (!pdbId.equals("1A0S")) continue; // C3
//			if (!pdbId.equals("1A5K")) continue; // C3
//			if (!pdbId.equals("1B44")) continue; // C5
//			if (!pdbId.equals("2Y9J")) continue; // C24 ?
//			if (!pdbId.equals("4HHB")) continue; // C2/D2
//			if (!pdbId.equals("1A95")) continue; // D2
//			if (!pdbId.equals("1A3G")) continue; // D3 good example for showing sym. axes
//			if (!pdbId.equals("1Q2V")) continue; // D8 // problem: C2 axes are too long ??
//			if (!pdbId.equals("2WCD")) continue; // D12
//			if (!pdbId.equals("1AEW")) continue; // O // interesting case, look for Cd along symmetry axes
//			if (!pdbId.equals("1A34")) continue; // I // dito 2x principal axes????
//			if (!pdbId.equals("1M4X")) continue; // I
//			if (!pdbId.equals("1COA")) continue; // D6
//			if (!pdbId.equals("1A5K")) continue; // dot= -0.9
//			if (!pdbId.equals("1M5Q")) continue; 
//			if (!pdbId.equals("3LSV")) continue;
//			if (!pdbId.equals("2BG9")) continue; // acetylcholin receptor, 2 alpha, 1 beta, 1 delta, 1 gamma
			// 1B4N, 1AVO, 1A5K, 1BZ5, 1A8S, 1BHC(D5), 1M5Q(1), 1KN0(C11) // good example
//			if (!pdbId.equals("1A6D")) continue; // trace <= 0 
//			if (!pdbId.equals("3KO1")) continue; // good
//			if (!pdbId.equals("2AO9")) continue; // good
//			if (!pdbId.equals("1NF4")) continue; // ???? octahedral//
//			if (!pdbId.equals("1H2I")) continue; // good

//			if (!pdbId.equals("1M5Q")) continue; // good example

			System.out.println("------------- " + pdbId  + "-------------");

		//	int bioAssemblyCount = entry.getBioAssemblyCount();
			int bioAssemblyCount = 1; // do only first bioassembly
			System.out.println("Bioassemblies: " + bioAssemblyCount);
			int n = Math.max(bioAssemblyCount, 1);
			for (int i = 0; i < n; i++) {		
				System.out.println("Bioassembly: " + i);
				Structure structure = null;
				int bioassemblyId = 0;
				try {
					if (bioAssemblyCount == 0) {
						structure = cache.getStructure(pdbId);
					} else {
						structure = cache.getBiologicalAssembly(pdbId, i+1, true);
						bioassemblyId = i+1;
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

				long tc1 = System.nanoTime(); 	

				QuatSymmetryParameters params = new QuatSymmetryParameters();
				params.setMinimumSequenceLength(MIN_SEQUENCE_LENGTH);
				params.setSequenceIdentityThreshold(SEQUENCE_IDENTITY_THRESHOLD);
				params.setAlignmentFractionThreshold(ALIGNMENT_FRACTION_THRESHOLD);
				params.setRmsdThreshold(RMSD_THRESHOLD);
				
				FindQuarternarySymmetry finder = new FindQuarternarySymmetry(structure, params);
				
				if (finder.getChainCount() == 0) {
					continue;
				}
				RotationGroup rotationGroup = finder.getRotationGroup();	
				String pointGroup = rotationGroup.getPointGroup();
				System.out.println("Point group: " + pointGroup);
				
				String formula = finder.getCompositionFormula();
				System.out.println("Formula: " + formula);
	

				// get metrics
				int caCount = finder.getSubunits().getCalphaCount();
				boolean groupComplete = rotationGroup.isComplete();
				float rmsd = (float) rotationGroup.getAverageSubunitRmsd();
				float rmsdT = (float) rotationGroup.getAverageTraceRmsd();
		
				int order = rotationGroup.getOrder();
				
				Subunits subunits = finder.getSubunits();
				List<String> chainIds = finder.getChainIds();
				int chainCount = subunits.getCenters().size();
				AxisTransformation at = new AxisTransformation(subunits, rotationGroup);
				float rx = (float) at.getDimension().x;
				float ry = (float) at.getDimension().y;
				float rz = (float) at.getDimension().z;
				int bin2 = Math.round(rz)/2 + 1000*(Math.round(ry)/2) + 1000000*(Math.round(rx)/2);
				Matrix4d matrix = at.getTransformation();
				JmolSymmetryScriptGenerator g = JmolSymmetryScriptGenerator.getInstance(at, rotationGroup);
                String jmolTransform = g.setOrientation(0);
				String jmolAxes = "";
				
				// determine overall symmetry
			
				MomentsOfInertia m = subunits.getMomentsOfInertia();
				MomentsOfInertia.SymmetryClass symmetryClass = m.getSymmetryClass(0.2);
				double asymmetryCoefficient = m.getAsymmetryParameter(0.001);
			

	//			if (true == false) {
				// create list of representative chains
	//			List<String> chainIds = finder.getChainIds();
				ProteinComplexSignature s100 = new ProteinComplexSignature(pdbId, chainIds, reader100);
				String signature100 = s100.getComplexSignature();
				String stoich100 = s100.getComplexStoichiometry();
				int types100 = s100.getSubunitTypeCount();
				ProteinComplexSignature s95 = new ProteinComplexSignature(pdbId, chainIds, reader95);
				String signature95 = s95.getComplexSignature();
				String stoich95 = s95.getComplexStoichiometry();
				int types90 = s95.getSubunitTypeCount();
				System.out.println("sign100: " + signature100 + ":" + stoich100);
				ProteinComplexSignature s70 = new ProteinComplexSignature(pdbId, chainIds, reader70);
				String signature70 = s70.getComplexSignature();
				String stoich70 = s70.getComplexStoichiometry();
				int types70 = s70.getSubunitTypeCount();
				System.out.println("sign70: " + signature70 + ":" + stoich70);
				ProteinComplexSignature s40 = new ProteinComplexSignature(pdbId, chainIds, reader40);
				String signature40 = s40.getComplexSignature();
				String stoich40 = s40.getComplexStoichiometry();
				int types40 = s40.getSubunitTypeCount();
				System.out.println("sign40: " + signature40 + ":" + stoich40);
				
				// TODO use chain signatures to label interactions of ligands
				long tc2 = System.nanoTime();
				long time = (tc2 - tc1)/1000000;
				symTime += time;


				
				// write .csv summary file
				if (chainCount >= MIN_CHAINS && groupComplete) {			
					out.print(pdbId + "," + bioassemblyId + "," + formula + "," + signature100 + "," + stoich100 + "," + types100 + "," + signature95 + "," + stoich95 + "," + types90 + "," + signature70 + "," + stoich70  + "," + types70 + "," + signature40 + "," + stoich40 + "," + types40 + "," + pointGroup + "," +
							order + "," + caCount + "," + chainCount + "," + rmsd + "," + rmsdT + "," + symmetryClass + "," + asymmetryCoefficient + "," +
						    time + "," + rx + ","  + ry + "," + rz + "," + bin2 + "," + jmolTransform + "," + jmolAxes + ",");
							out.println();
							out.flush();
							if (chainCount > 1) {
								multimer++;
							}
				} else if (chainCount >= MIN_CHAINS) {
					out1.print(pdbId + "," + bioassemblyId + "," + formula + "," + signature100 + "," + stoich100 + "," + types100 + "," + signature95 + "," + stoich95 + "," + types90 + "," + signature70 + "," + stoich70  + "," + types70 + "," + signature40 + "," + stoich40 + "," + types40 + "," + pointGroup + "," + groupComplete + "," +
							order + "," + caCount + "," + chainCount  + "," + rmsd + "," + rmsdT + "," + symmetryClass + "," + asymmetryCoefficient + "," +
						    time + "," + rx + ","  + ry + "," + rz + "," + bin2 + "," + jmolTransform + "," + jmolAxes + ",");
					out1.println();
					out1.flush();
					excluded++;
				}
	//			}

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
}
