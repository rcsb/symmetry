package org.biojava3.structure.quaternary.misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;

public class FindPdbRepresentatives {
	private static final String SERVICELOCATION="http://www.rcsb.org/pdb/rest/postBLAST";
	private static double eCutoffValue = 0.00001;
	
	private Structure structure = null;
	
	public FindPdbRepresentatives(Structure structure) {
		this.structure = structure;
	}
	
	public List<PdbBlastHit> findBestBlastHits() throws Exception {
		List<PdbBlastHit> representatives = new ArrayList<PdbBlastHit>();
		
		for (String sequence: getUniqueSequences()) {
			List<PdbBlastHit> hits = runBlast(sequence);	
			PdbBlastXMLParser.sortBlastHits(hits);
			
//			System.out.println("Sorted hits");
//			for (PdbBlastHit hit: hits) {
//				System.out.println(hit);
//			}
			if (! hits.isEmpty()) {
		  	    representatives.add(hits.get(0));
			}
		}
		return representatives;
	}
	
	private List<String> getUniqueSequences() {
		List<String> sequences = new ArrayList<String>();
		for (Chain c: structure.getChains()) {
			String seq = c.getSeqResSequence();
			// TODO only add amino acid seq.
			if (! sequences.contains(seq)) {
				sequences.add(seq);
			}
		}
		return sequences;
	}
	
	private List<PdbBlastHit> runBlast(String param1) {
		String param2 = "eCutOff=" + eCutoffValue;     
		String param3 = "matrix=BLOSUM62"; 
		String param4 = "outputFormat=XML";  // HTML or XML. If not specified, default to plain text 

		try {
			// Send the request 
			URL url = new URL(SERVICELOCATION);
			URLConnection conn = url.openConnection();
			conn.setDoOutput(true); 
			BufferedWriter out = new BufferedWriter(new OutputStreamWriter(conn.getOutputStream()));

			// Write parameters 
			out.write("sequence=");
			out.write(param1);
			out.write("&");
			out.write(param2);
			out.write("&");
			out.write(param3);
			out.write("&");
			out.write(param4);
			out.flush();
			out.close();

			 // Get the response
//	         StringBuffer answer = new StringBuffer();
//	         BufferedReader in = new BufferedReader( new InputStreamReader( conn.getInputStream()) );
//	         String line;
//	         while ( (line = in.readLine()) != null ) {
//	            System.out.println(line);
//	         }
//	         in.close();
//	            
//	         // Output the response
//	         System.out.println(answer.toString());
			
			// Get the response
			PdbBlastXMLParser parser = new PdbBlastXMLParser(conn.getInputStream());
			return parser.parse(eCutoffValue);
		}
		catch (Exception ex) {
			ex.printStackTrace();
		}

		return Collections.emptyList();
	}
}
