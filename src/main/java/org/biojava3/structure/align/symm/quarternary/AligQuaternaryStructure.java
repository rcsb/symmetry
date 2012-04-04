package org.biojava3.structure.align.symm.quarternary;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;

public class AligQuaternaryStructure {
	private static String PDB_PATH = "C:/PDB/";
	private static final int MIN_SEQUENCE_LENGTH = 24;
	private Map<Integer,Integer> clusterMap = new LinkedHashMap<Integer,Integer>();
    private ClusterAlignment calign = new ClusterAlignment();
	
	public double[] align(String pdbId1, String pdbId2) {
		double[] bestRmsd = {Double.MAX_VALUE, Double.MAX_VALUE};
		Structure s1 = getStructure(pdbId1);
		if (s1 == null) {
			return bestRmsd;
		}
		GlobalSequenceGrouper grouper1 = new GlobalSequenceGrouper(s1, MIN_SEQUENCE_LENGTH);

		Structure s2 = getStructure(pdbId2);
		if (s2 == null) {
			return bestRmsd;
		}
		GlobalSequenceGrouper grouper2 = new GlobalSequenceGrouper(s2, MIN_SEQUENCE_LENGTH);

		ClusterAlignment alignment = alignUniqueSequences(grouper1, grouper2);

		Matrix4d bestTransformation = new Matrix4d();
		SubunitMapper mapper = new SubunitMapper(s1, s2);
		List<List<Integer>> mappings = mapper.getMappings();
		//		List<List<Integer>> mappings = mapper.getAllMappings();

		for (List<Integer> mapping: mappings) {
			Matrix4d transformation = new Matrix4d();
			double[] rmsds = quaternaryStructureAlignment(mapping, alignment,grouper1, grouper2, transformation);
			System.out.println("Mapping: " + mapping + " RMSD: " + rmsds[0]);
			if (rmsds[0] < bestRmsd[0]) {
				bestRmsd[0] = rmsds[0];
				bestRmsd[1] = rmsds[1];
				bestTransformation.set(transformation);
			}
		}

		System.out.println("RESULT: rmsd: " + bestRmsd);
		QuatSymmetryWriter writer = new QuatSymmetryWriter(s1);
		String fileName = PDB_PATH + "transformed/" + pdbId1 + ".pdb";

		writer.writeTransformedStructure(bestTransformation, fileName);
		//				writer.writeTransformedStructureNotCentered(m, fileName);
//		//				
		fileName = PDB_PATH + "transformed/" + pdbId2 + ".pdb";
		Matrix4d identity = new Matrix4d();
		identity.setIdentity();
		writer = new QuatSymmetryWriter(s2);
		writer.writeTransformedStructure(identity, fileName);


		return bestRmsd;
	}
	
	public double[] align(Structure s1, Structure s2) {
		double[] bestRmsd = {Double.MAX_VALUE, Double.MAX_VALUE};
		if (s1 == null) {
			return bestRmsd;
		}
		GlobalSequenceGrouper grouper1 = new GlobalSequenceGrouper(s1, MIN_SEQUENCE_LENGTH);

		if (s2 == null) {
			return bestRmsd;
		}
		GlobalSequenceGrouper grouper2 = new GlobalSequenceGrouper(s2, MIN_SEQUENCE_LENGTH);

		ClusterAlignment alignment = alignUniqueSequences(grouper1, grouper2);

		Matrix4d bestTransformation = new Matrix4d();
		SubunitMapper mapper = new SubunitMapper(s1, s2);
		List<List<Integer>> mappings = mapper.getMappings();
		//		List<List<Integer>> mappings = mapper.getAllMappings();

		for (List<Integer> mapping: mappings) {
			Matrix4d transformation = new Matrix4d();
			double[] rmsds = quaternaryStructureAlignment(mapping, alignment,grouper1, grouper2, transformation);
			System.out.println("Mapping: " + mapping + " RMSD: " + rmsds[0]);
			if (rmsds[0] < bestRmsd[0]) {
				bestRmsd[0] = rmsds[0];
				bestRmsd[1] = rmsds[1];
				bestTransformation.set(transformation);
			}
		}

		System.out.println("RESULT: rmsd: " + bestRmsd);
//		QuatSymmetryWriter writer = new QuatSymmetryWriter(s1);
//		String fileName = PDB_PATH + "transformed/" + pdbId1 + ".pdb";
//
//		writer.writeTransformedStructure(bestTransformation, fileName);
//		//				writer.writeTransformedStructureNotCentered(m, fileName);
//		//				
//		fileName = PDB_PATH + "transformed/" + pdbId2 + ".pdb";
//		Matrix4d identity = new Matrix4d();
//		identity.setIdentity();
//		writer = new QuatSymmetryWriter(s2);
//		writer.writeTransformedStructure(identity, fileName);


		return bestRmsd;
	}
	
	private double[] quaternaryStructureAlignment(List<Integer> mapping, ClusterAlignment alignment, GlobalSequenceGrouper grouper1, GlobalSequenceGrouper grouper2, Matrix4d transformation) {
		List<Atom[]> cas1 = grouper1.getCalphaTraces();
		List<Atom[]> cas2 = grouper2.getCalphaTraces();
		List<Atom[]> cas1a = new ArrayList<Atom[]>();
		List<Atom[]> cas2a = new ArrayList<Atom[]>();
		double[] rmsds = new double[2];
		
		double trmsd = 0;
		double tsum = 0;
		
		for (int i = 0; i < cas1.size(); i++) {
			int subunitId = mapping.get(i);
//			System.out.println("Alignment len1: " + alignment.getAlignment1(i).size());
			Atom[] ca1 =  createCalphaList(cas1.get(i), alignment.getAlignment1(i));
			cas1a.add(ca1);

//			System.out.println("Alignment len2: " + alignment.getAlignment2(i).size());
			Atom[] ca2 =  createCalphaList(cas2.get(subunitId), alignment.getAlignment2(i));
			if (ca1 == null || ca2 == null || ca1.length != ca2.length) {
				System.err.println("quaternaryStructureAlignment: ERROR: sequence length mismatch");
				rmsds[0] = 99999;
				rmsds[1] = 99999;
				return rmsds;
			}
			cas2a.add(ca2);	
			
			Point3d[] tpoints1 = getPoints(ca1);
			Point3d[] tpoints2 = getPoints(ca2);
			Matrix4d m = SuperPosition.superposeWithTranslation(tpoints1, tpoints2);
			double sum = SuperPosition.rmsd(tpoints1, tpoints2);
			// "undo" RMSD to get raw sum of distance squares
			sum *= sum;
			sum *= tpoints1.length;
			tsum += sum;
		}
		Point3d[] points1 = getPoints(cas1a);
		Point3d[] points2 = getPoints(cas2a);
		Matrix4d m = SuperPosition.superposeWithTranslation(points1, points2);
		double qrmsd = SuperPosition.rmsd(points1, points2);
//		System.out.println("Quaternary structure rmsd: " + rmsd);
		transformation.set(m);
		trmsd = Math.sqrt(tsum/points1.length);
		
		rmsds[0] = qrmsd;
		rmsds[1] = trmsd;
		return rmsds;		
	}
	
	/**
	 * Converts a list of C-alpha traces to an array of C-alpha coordinates
	 * @param traces
	 * @return
	 */
	private Point3d[] getPoints(List<Atom[]> traces) {
		int len = 0;
		for (Atom[] atoms: traces) {
			len += atoms.length;
		}
		Point3d[] points = new Point3d[len];
		int i = 0;
		for (Atom[] atoms: traces) {
			for (Atom a: atoms) {
				points[i] = new Point3d(a.getCoords());
				i++;
			}
		}
		return points;
	}
	
	private Point3d[] getPoints(Atom[] trace) {
		Point3d[] points = new Point3d[trace.length];
		for (int i = 0; i < trace.length; i++) {
			points[i] = new Point3d(trace[i].getCoords());
		}
		return points;
	}
	
	private Atom[] createCalphaList(Atom[] caSeq, List<Integer> alignment) {
		Atom[] subset = new Atom[alignment.size()];
//		System.out.println("createCalphaList: len/alignment: " + caSeq.length + ": "+ alignment);
		if (alignment.size() > caSeq.length) {
			System.out.println("createCalphaList: ERROR: size mismatch");
			return null;
		}
		for (int i = 0; i < alignment.size(); i++) {
			if (alignment.get(i) >= caSeq.length) {
				return null;
			}
			subset[i] = caSeq[alignment.get(i)];
		}
		return subset;
	}
	
	private ClusterAlignment alignUniqueSequences(GlobalSequenceGrouper grouper1, GlobalSequenceGrouper grouper2) {
		List<List<Integer>> clusters1 = grouper1.getSequenceCluster100();
	    List<Atom[]> ca1 = grouper1.getCalphaTraces();
  //      System.out.println("Cluster1: " + clusters1.size());

        List<List<Integer>> clusters2 = grouper2.getSequenceCluster100();
        List<Atom[]> ca2 = grouper2.getCalphaTraces();
 //       System.out.println("Cluster2: " + clusters2.size());
        clusterMap = new LinkedHashMap<Integer,Integer>();
  
        calign = new ClusterAlignment();
        
        boolean first = true;
        
        for (List<Integer> c1: clusters1) {
        	// get the first C alpha array from a cluster, as the representative of this cluster
        	int representative1 = c1.get(0); 

        	Atom[] ca1Seq = ca1.get(representative1);
        	for (List<Integer> c2: clusters2) {
        		// can't match cluster of different sizes
        		if (c1.size() != c2.size()) {
        			if (first) {
        				System.out.println("WARNING: cluster 1 and 2 don't match");
        			}
        			continue;
        		}
        		
        		int representative2 = c2.get(0); 
        		if (clusterMap.containsValue(representative2)) {
        			continue;
        		}
            	System.out.print("Comparing: " + representative1);
            	System.out.println(" - " + representative2);
            	Atom[] ca2Seq = ca2.get(representative2);
            	System.out.println("Ca length: " + ca1Seq.length + " - " + ca2Seq.length);
        		AFPChain afp = alignPair(ca1Seq, ca2Seq);
        		if (afp == null) {
        			System.out.println("AFPChain is null");
        			continue;
        		}
        		double identity = afp.getIdentity();
        		
        		if (first && identity <= 0.95) {
    				System.out.println("WARNING: cluster 1 and 2 have insufficient seq. id: " + identity);
    				continue;
    			}
        		first = false;
        		// store alignment information for two matching clusters
        		// TODO should only keep best match, id may be as low at 0.9
        		if (identity > 0.95) {
        			System.out.println("Seq. identity: " + afp.getIdentity());
        			System.out.println("Rmsd: " + afp.getChainRmsd());
        			int[][][] alig = afp.getOptAln();
 //       			System.out.println("Align1: " + Arrays.toString(alig[0][0]));
 //       			System.out.println("Align2: " + Arrays.toString(alig[0][1]));
        			List<Integer> alignment1 = new ArrayList<Integer>();
        			if (alig == null) {
        				continue;
        			}
        			for (Integer a1: alig[0][0]) {
        				alignment1.add(a1);
        			}
        			List<Integer> alignment2 = new ArrayList<Integer>();
        			for (Integer a2: alig[0][1]) {
        				alignment2.add(a2);
        			}
        			if (alignment1.size() != alignment2.size()) {
        				System.out.println("alignUniqueSequences: ERROR: alignment length mismatch");
        				continue;
        			}
        			calign.add(c1, c2, alignment1, alignment2, identity);
        			clusterMap.put(representative1,representative2);
        			break;
        		}
        		
        		
        	}
        }
        return calign;
	}
	
	private AFPChain alignPair(Atom[] ca1Seq, Atom[] ca2Seq) {
		SmithWaterman3Daligner aligner = new SmithWaterman3Daligner();
		AFPChain afp = null;
		try {
			afp = aligner.align(ca1Seq, ca2Seq);
		} catch (StructureException e) {
			e.printStackTrace();
			return afp;
		} 
		return afp;
	}
	
	private Structure getStructure(String pdbId) {
		AtomCache cache = new AtomCache();
		cache.setAutoFetch(true);
		cache.setPath(PDB_PATH);	
		FileParsingParameters p = new FileParsingParameters();
		p.setStoreEmptySeqRes(true);
		//p.setParseCAOnly(true);
		//p.setMaxAtoms(50000000);
		
		p.setAtomCaThreshold(Integer.MAX_VALUE);
	//	p.setAcceptedAtomNames(new String[]{" CA ", " CB "});
		cache.setFileParsingParams(p);
			
		PrintWriter error = null;
		try {
			error = new PrintWriter(new FileWriter(PDB_PATH + "error.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		Structure structure = null;
		try {
			structure = cache.getBiologicalAssembly(pdbId, 1, true);
			//				structure = cache. getStructure(pdbId);

		} catch (StructureException e) {
			e.printStackTrace();
			error.println(pdbId + "------------------------------------");
			error.println(e.getMessage());
			error.flush();
		} catch (IOException e) {
			e.printStackTrace();
			error.println(pdbId + "------------------------------------");
			error.println(e.getMessage());
			error.flush();
		}
		return structure;
	}
}
