/**
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on Mar 7, 2013
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package demo;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.symmetry.analysis.CalcBioAssemblySymmetry;
import org.biojava.bio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.bio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava3.structure.StructureIO;

public class TestLoadingManyStoichs {
	public static final String SERVICELOCATION="http://www.rcsb.org/pdb/rest";
	
	
	public static final String testPointGroup="D2";

	public static void main(String[] args){
		System.out.println("loading pointGroups with " + testPointGroup);

		TestLoadingManyStoichs me = new TestLoadingManyStoichs();

		List<String> pdbIDs = me.getPDBIDs();

		for (String pdbID : pdbIDs){
			calculateStoichiometry(pdbID);
		}
	}


	private static void calculateStoichiometry(String pdbID) {
		if (pdbID == null || pdbID.length() != 4)
			return;

		System.out.println("calculating symmetry for " + pdbID);

		runPDB(pdbID);

	}

	public static void runPDB(String pdbID){

		pdbID = pdbID.toLowerCase();

		int  biolAssemblyNr =1;

		Structure s;
		try {

			//			
			AtomCache cache = new AtomCache();
			FileParsingParameters params = cache.getFileParsingParams();
			params.setAlignSeqRes(true);
			params.setParseCAOnly(false);

			StructureIO.setAtomCache(cache);

			s = StructureIO.getBiologicalAssembly(pdbID, biolAssemblyNr);

			// Alternative access to structure:			
			//
			//s = readStructure(pdbID, biolAssemblyNr);

			//System.out.println("MODELS:" + s.nrModels());
			analyzeSymmetry(s,pdbID, 0, 0.30);
			String  symmetry = analyzeSymmetry(s,pdbID, biolAssemblyNr, 0.30);

			System.out.println(" *** " + pdbID + " " + biolAssemblyNr + " " + symmetry);			
			if ( ! symmetry.equals(testPointGroup)) {
				System.exit(0);
			}

			symmetry = analyzeSymmetry(s,pdbID, biolAssemblyNr, 0.95);

			System.out.println(" *** " + pdbID + " " + biolAssemblyNr + " " + symmetry);	
			if ( ! symmetry.equals(testPointGroup)) {
				System.exit(0);
			}



		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}

	private static String analyzeSymmetry(Structure s,String pdbID, int biolAssemblyNr, double threshold) {
        QuatSymmetryParameters parameters = new QuatSymmetryParameters();
	    parameters.setVerbose(false);

		CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry(s, parameters);
		QuatSymmetryDetector detector = calc.orient();
		boolean hasProtein = detector.hasProteinSubunits();

		String symmetry = null;
		if ( hasProtein) {
			symmetry = calc.getRotationGroup().getPointGroup().toString();
		}

		return symmetry;
	}

	public List<String> getPDBIDs(){

		String xml =
				"<orgPdbCompositeQuery version=\"1.0\">"+
						" <queryRefinement>"+
						"  <queryRefinementLevel>0</queryRefinementLevel>"+
						"		 <orgPdbQuery>"+
						"		    <version>head</version>"+
						"		    <queryType>org.pdb.query.simple.PointGroupQuery</queryType>"+
						"		    <description>Finds PDB entries based on symmetry: Point Group is T and Sequence ID is 0.95</description>"+
						"		    <queryId>E316FC2F</queryId>"+
						"		    <resultCount>257</resultCount>"+
						"		    <runtimeStart>2013-03-07T22:03:23Z</runtimeStart>"+
						"		    <runtimeMilliseconds>100</runtimeMilliseconds>"+
						"		    <pointGroup>"+ testPointGroup +"</pointGroup>"+
						"		    <sequenceID>0.95</sequenceID>"+
						"		  </orgPdbQuery>"+
						" </queryRefinement>"+
						"</orgPdbCompositeQuery>" ;


		List<String> pdbIds = new ArrayList<String>();
		try {
			pdbIds = postQuery(xml);

		} catch (Exception e){
			e.printStackTrace();
		}

		return pdbIds;
	}


	/** post am XML query (PDB XML query format)  to the RESTful RCSB web service
	 * 
	 * @param xml
	 * @return a list of PDB ids.
	 */
	public List<String> postQuery(String xml) 
			throws IOException{


		URL u = new URL(SERVICELOCATION + "/search");

		String encodedXML = URLEncoder.encode(xml,"UTF-8");

		InputStream in =  doPOST(u,encodedXML);

		List<String> pdbIds = new ArrayList<String>();

		BufferedReader rd = new BufferedReader(new InputStreamReader(in));
		String line;
		while ((line = rd.readLine()) != null) {
			pdbIds.add(line);

		}      
		rd.close();

		return pdbIds;



	}

	/** do a POST to a URL and return the response stream for further processing elsewhere.
	 * 
	 * 
	 * @param url
	 * @return
	 * @throws IOException
	 */
	public static InputStream doPOST(URL url, String data)
			throws IOException 
	{

		// Send data

		URLConnection conn = url.openConnection();
		conn.setDoOutput(true);

		OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
		wr.write(data);
		wr.flush();

		// Get the response
		return conn.getInputStream();

	}
}
