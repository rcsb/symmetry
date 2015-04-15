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
 * Created on Oct 5, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.List;

public class FileUtils {
	public static String readFileAsString(String filePath)
	throws java.io.IOException{
		StringBuffer fileData = new StringBuffer(1000);
		BufferedReader reader = new BufferedReader(
				new FileReader(filePath));
		char[] buf = new char[1024];
		int numRead=0;
		while((numRead=reader.read(buf)) != -1){
			String readData = String.valueOf(buf, 0, numRead);
			fileData.append(readData);
			buf = new char[1024];
		}
		reader.close();
		return fileData.toString();
	}
	
	public static void writeStringToFile(String content, String filePath){
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter( new FileWriter( filePath));
			writer.write( content);
		}
		catch ( IOException e)
		{
			e.printStackTrace();
			
		}
		finally
		{
			try
			{
				if ( writer != null)
					writer.close( );
			}
			catch ( IOException e)
			{
			}
		}

		
	}
	
	public static void copy(File src, File dst) throws IOException {

		InputStream in = new FileInputStream(src);
		OutputStream out = new FileOutputStream(dst);

		// Transfer bytes from in to out
		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0) {
			out.write(buf, 0, len);
		}
		in.close();
		out.close();
	}

	/**
	 * Saves a graph into a csv file in the format of tuples (vertex,edge) for every edge in the graph.
	 * The graph has to be in the form of an adjacency list.
	 */
	 public static void saveGraph(List<List<Integer>> graph, String sFileName) throws IOException {
		
	    FileWriter writer = new FileWriter(sFileName);
	    writer.append("Vertex,Edge\n");
	    for (Integer i=0; i<graph.size(); i++){
	    	for (int j=0; j<graph.get(i).size(); j++){
	    		
	    		writer.append(i.toString());
	    		writer.append(',');
	    		writer.append(graph.get(i).get(j).toString());
	    		writer.append('\n');
	    	}
	    }
	    
	    writer.flush();
	    writer.close();
	 }
	 
	/**
	 *  Save a multiple sequence alignment result as a fasta file.
	 */
	public static void saveSeqAln(String[] alnSequences, String filePath) throws IOException {
		
		FileWriter writer = new FileWriter(filePath);
		
		for (int su=0; su<alnSequences.length; su++){
			writer.append(">Subunit_"+(su+1)+"\n");
			writer.append(alnSequences[su]+"\n");
		}
		
		writer.flush();
	    writer.close();
	}

}
