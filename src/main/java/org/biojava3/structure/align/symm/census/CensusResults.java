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
 * Created on Sep 30, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava3.structure.align.symm.census;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.utils.FileUtils;




@XmlRootElement(name = "CensusResults", namespace ="http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class CensusResults implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -5517546595033480440L;
	static JAXBContext jaxbContext;
	static {
		try {
			jaxbContext= JAXBContext.newInstance(CensusResults.class);
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	private List<CensusResult> data = new ArrayList<CensusResult>();
	
	
	
	public List<CensusResult> getData() {
		return data;
	}

	public void setData(List<CensusResult> data) {
		this.data = data;
	}

	public String toHTML(){
		StringWriter w = new StringWriter();
		
		w.append(ResultConverter.getHTMLHeader());
		for ( CensusResult d : data){
			w.append(ResultConverter.toHTML(d));
		}
		
		w.append("</tbody></table>");
		return w.toString();
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

	public static CensusResults fromXML(String xml){

		CensusResults job = null;

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();

			ByteArrayInputStream bais = new ByteArrayInputStream(xml.getBytes());

			job = (CensusResults) un.unmarshal(bais);

		} catch (Exception e){
			e.printStackTrace();
		}

		return job;
	}
	
	public static CensusResults getExistingResults() throws IOException{
		AtomCache cache  = new AtomCache();

		String r = cache.getPath() + File.separator + "scopCensus.xml";
		File f = new File(r); 
		if (f.exists()) {
			
			String xml = FileUtils.readFileAsString(r);
			CensusResults data = CensusResults.fromXML(xml);
			
			System.out.println("read " + data.getData().size() + " results from disk...");
			return data;
		}
		return null;
	}
	

}
