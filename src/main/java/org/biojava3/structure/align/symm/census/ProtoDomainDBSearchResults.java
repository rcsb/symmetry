package org.biojava3.structure.align.symm.census;

import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

import org.biojava3.structure.utils.FileUtils;



/**
 * @deprecated Refer to {@code symm.census2} instead
 */
@Deprecated
@XmlRootElement(name = "ProtoDomainDBSearchResults", namespace ="http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class ProtoDomainDBSearchResults implements Serializable {


	/**
	 * 
	 */
	private static final long serialVersionUID = -5517546595033480440L;
	static JAXBContext jaxbContext;
	static {
		try {
			jaxbContext= JAXBContext.newInstance(ProtoDomainDBSearchResults.class);
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	private List<ProtoDomainDBSearchResult> dbSearchResults;
	String protoDomain;

	public static ProtoDomainDBSearchResults fromXML(String xml) {


		ProtoDomainDBSearchResults job = null;

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();

			ByteArrayInputStream bais = new ByteArrayInputStream(xml.getBytes());

			job = (ProtoDomainDBSearchResults) un.unmarshal(bais);

		} catch (Exception e){
			e.printStackTrace();
		}

		return job;
	}
	public  String toXML(){

		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		PrintStream ps = new PrintStream(baos);

		try {

			Marshaller m = jaxbContext.createMarshaller();

			m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

			m.marshal( this, ps);
			

		} catch (Exception e){
			e.printStackTrace();
		}

		return baos.toString();

	}
	public List<ProtoDomainDBSearchResult> getDbSearchResults() {
		return dbSearchResults;
	}
	public void setDbSearchResults(List<ProtoDomainDBSearchResult> dbSearchResults) {
		this.dbSearchResults = dbSearchResults;
	}
	
	
	public static ProtoDomainDBSearchResults fromFile(File resultsFile) throws IOException{
		String dbResultsXML = FileUtils.readFileAsString(resultsFile.getAbsolutePath());

		ProtoDomainDBSearchResults dbResults = ProtoDomainDBSearchResults.fromXML(dbResultsXML);
		
		return dbResults;

	}
	
	public void writeToFile(File resultsFile) {
		
		String xml = toXML();
		
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter( new FileWriter( resultsFile));
			writer.write( xml);
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

		System.out.println(resultsFile.getAbsolutePath() + " wrote " + dbSearchResults.size() + " results to disk...");
		
	}
	public String getProtoDomain() {
		return protoDomain;
	}
	public void setProtoDomain(String protoDomain) {
		this.protoDomain = protoDomain;
	}



}
