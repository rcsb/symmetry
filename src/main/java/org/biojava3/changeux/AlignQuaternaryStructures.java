package org.biojava3.changeux;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.SortedSet;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava3.core.util.InputStreamProvider;

import org.rcsb.fatcat.server.PdbChainKey;
import org.rcsb.fatcat.server.dao.PdbDAO;
import org.rcsb.fatcat.server.util.SetupJNDIDataSource;

public class AlignQuaternaryStructures {

	protected static AtomCache cache = new AtomCache("/Users/ap3/WORK/PDB/",true);

	private static final String separator = System.getProperty("line.separator");
	
	public static void main(String[] args){



		int clusterCutoff = 95;
		//SortedSet<PdbChainKey>  representatives = getSequenceClusters(clusterCutoff);

		//for ( PdbChainKey repre: representatives) {

		PdbChainKey repre =  PdbChainKey.fromName("3KFK.A");
		//PdbChainKey repre =  PdbChainKey.fromName("4HHB.A");

		try {
			alignCluster(repre,clusterCutoff);
		} catch (Exception e){
			e.printStackTrace();
		}
		//}

	}

	private static void alignCluster(PdbChainKey repre,int clusterCutoff) 
	throws IOException, StructureException {



		PdbDAO dao = new PdbDAO();
		int clusterId = dao.getClusterNumber(clusterCutoff,repre);
		SortedSet<PdbChainKey> members = dao.getClusterMembers(clusterCutoff,clusterId);

		System.out.println("### cluster " + clusterCutoff + " has " + members.size() + " members");


		StructureAlignment algo = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);

		//CeParameters params = (CeParameters)algo.getParameters();
		//params.setMaxGapSize(-1);

		File  logFile = new File("/tmp/repre_" + repre.getPdbId()+".txt");
		
		for (PdbChainKey member : members){
			String xml = null;
			String outputFileName = "/Users/ap3/WORK/PDB/bio/"+ repre.getPdbId() + "_" + member.getPdbId()+".xml.gz";
			File f = new File(outputFileName);
			if ( f.exists()){
				System.out.println("loading results from previous calculation for " + f);
				xml = showPrecalcResult(repre, member,f, logFile);
			} else {

				
				try {
					xml = align(repre, member,algo,logFile);
					writeXML2File(f, xml);
				} catch (Exception e){
					//e.printStackTrace();
					System.err.println(e.getMessage());
				}
			}
		}


	}
	
	private static void writeToResultsFile(File f, String msg){
		try{
			// Create file 
			FileWriter fstream = new FileWriter(f,true);
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(msg);
			out.write(separator);
			//Close the output stream
			out.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}

	private static String showPrecalcResult(PdbChainKey repre,
			PdbChainKey member, File target, File logFile) throws StructureException {
		
		String xml = null;
		try {
			xml = getXMLFromFile(repre, member, target);
			
			Structure s1 = cache.getBiologicalUnit(repre.getPdbId());
			Structure s2 = cache.getBiologicalUnit(member.getPdbId());

			s1.getPDBHeader().setTitle("Biological Unit of " + repre.getPdbId());
			s2.getPDBHeader().setTitle("Biological Unit of " + member.getPdbId());
			s1.setPDBCode(repre.getPdbId());
			s2.setPDBCode(member.getPdbId());
			Atom[] ca1 = StructureTools.getAtomCAArray(s1);
			Atom[] ca2 = StructureTools.getAtomCAArray(s2);
			
			AFPChain afpChain = AFPChainXMLParser.fromXML(xml, ca1,ca2);
			double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
			afpChain.setTMScore(tmScore);
			
			String msg = repre.toName() + "\t" + member.toName() + "\t" + 
			afpChain.getProbability() + "\t" + afpChain.getCoverage1() + "\t" + 
			afpChain.getCoverage2() +"\t" + tmScore+ "\t-1" ;
			writeToResultsFile(logFile, msg);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return xml;
	}

	private static String getXMLFromFile(PdbChainKey repre, PdbChainKey member, File target) throws IOException {
		  
		InputStreamProvider prov = new InputStreamProvider();
		InputStream in = prov.getInputStream(target);
		StringWriter sw = new StringWriter();

		BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(in));
		String line;
		while ((line = fileBuffer.readLine()) != null) {
			sw.append(line);
			sw.append(separator);
		}
		sw.flush();
		sw.close();
		return sw.toString();
	}

	public static String align(PdbChainKey repre,
			PdbChainKey member, StructureAlignment algo,
			File logFile
	) throws StructureException, IOException {
		Structure s1 = cache.getBiologicalUnit(repre.getPdbId());
		Structure s2 = cache.getBiologicalUnit(member.getPdbId());

		s1.getPDBHeader().setTitle("Biological Unit of " + repre.getPdbId());
		s2.getPDBHeader().setTitle("Biological Unit of " + member.getPdbId());
		s1.setPDBCode(repre.getPdbId());
		s2.setPDBCode(member.getPdbId());
		String[] atomNames = {StructureTools.caAtomName };
		Atom[] ca1 = StructureTools.getAtomArrayAllModels(s1,atomNames);
		Atom[] ca2 = StructureTools.getAtomArrayAllModels(s2,atomNames);

//		StructureAlignmentJmol jmol1 = new  StructureAlignmentJmol(null,null,null);
//		jmol1.setStructure(s1);
//
//		StructureAlignmentJmol jmol2 = new  StructureAlignmentJmol(null,null,null);
//		jmol2.setStructure(s2);


		String msg1 = "aligning " + s1.getPDBCode() + " size:" + ca1.length + 
				" atoms | " + s2.getPDBCode() + " size:" + ca2.length + " atoms";
		System.out.println(msg1);
		
				
		long startTime =System.currentTimeMillis();
		AFPChain afpChain = algo.align(ca1, ca2);
		
		double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
		afpChain.setTMScore(tmScore);
		
		long endTime = System.currentTimeMillis();
		
		String msg = repre.toName() + "\t" + member.toName() + "\t" + 
				afpChain.getProbability() + "\t" + afpChain.getCoverage1() + "\t" + 
				afpChain.getCoverage2() +"\t" + tmScore+ "\t" + (endTime - startTime)/1000 ;
		if ( logFile != null)
			writeToResultsFile(logFile, msg);
		else {
			System.out.println(msg);
		}
		
		String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);
		return xml;
	}

	private static void writeXML2File(File target, String xml) throws FileNotFoundException, IOException{
		GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(target));

		InputStream in = new ByteArrayInputStream(xml.getBytes("UTF-8"));    

		// Transfer bytes from the input file to the GZIP output stream
		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0) {
			out.write(buf, 0, len);
		}
		in.close();
		out.finish();
		out.close();
	}

	private static SortedSet<PdbChainKey> getSequenceClusters(int i) {

		PdbDAO dao = new PdbDAO();

		return dao.getRepresentatives(i);

	}


}
