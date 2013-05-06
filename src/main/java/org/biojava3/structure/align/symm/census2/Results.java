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
package org.biojava3.structure.align.symm.census2;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

/**
 * A census (list of census {@link Result Results}).
 * @author dmyerstu
 */
@XmlRootElement(name = "CensusResults", namespace = "http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class Results implements Serializable {

	private static final long serialVersionUID = -5517546595033480440L;
	private static JAXBContext jaxbContext;
	private List<Result> data = new ArrayList<Result>();

	private String timestamp;

	static {
		try {
			jaxbContext = JAXBContext.newInstance(Results.class);
		} catch (Exception e) {
			throw new RuntimeException(e); // fatal
		}
	}

	public static Results fromXML(File file) throws IOException {

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();
			FileInputStream fis = new FileInputStream(file);
			Results results = (Results) un.unmarshal(fis);

			// due to a side effect by JAXB
			List<Result> newData = new ArrayList<Result>(results.getData().size());
			
			for (Result result : results.getData()) {
				if (result != null) newData.add(result);
			}
			results.setData(newData);

			return results;

		} catch (JAXBException e) {
			throw new IOException(e);
		}

	}
	
	public static List<SimpleResult> convertResults(List<Result> data){
		List<SimpleResult> results = new ArrayList<SimpleResult>();
		
		for ( Result r : data){
			SimpleResult s = new SimpleResult();
			results.add(s);
			
			s.setAlignScore(r.getAlignment().getAlignScore());
			s.setClassification(r.getClassification());
			s.setDescription(r.getDescription());
			s.setGapLength(r.getAlignment().getGapLength());
			s.setIdentity(r.getAlignment().getIdentity());
			s.setIsSignificant(r.getIsSignificant());
			s.setOrder(r.getOrder());
			s.setProtodomain(r.getProtodomain());
			s.setRank(r.getRank());
			s.setScopId(r.getScopId());
			s.setSimilarity(r.getAlignment().getSimilarity());
			s.setSunId(r.getSunId());
			if ( r.getAxis() != null)
				s.setTheta(r.getAxis().getTheta());
			s.setTmScore(r.getAlignment().getTmScore());
			s.setzScore(r.getAlignment().getzScore());
			s.setRmsd(r.getAlignment().getRmsd());
			s.setAlignLength(r.getAlignment().getAlignLength());
			s.setLength1(r.getAlignment().getGapLength() + r.getAlignment().getAlignLength());
		}
		
		return results;
	}

	public static Results fromXML(File[] files) throws IOException {
		Results results = new Results();
		for (File file : files) {
			results.getData().addAll(fromXML(file).getData());
		}
		return results;
	}

	public static Results getExistingResults(File file) throws IOException {
		if (file.exists()) {
			return Results.fromXML(file);
		}
		return null;
	}

	public static Results getExistingResults(String file) throws IOException {
		return getExistingResults(new File(file));
	}

	public Results() {
		timestamp = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date());
	}

	public boolean add(Result e) {
		return data.add(e);
	}

	public boolean addAll(Collection<? extends Result> c) {
		return data.addAll(c);
	}

	public int size() {
		return data.size();
	}

	public boolean contains(Object o) {
		return data.contains(o);
	}

	public boolean containsAll(Collection<?> c) {
		return data.containsAll(c);
	}

	public List<Result> getData() {
		return data;
	}

	public String getTimestamp() {
		return timestamp;
	}

	public void setData(List<Result> data) {
		this.data = data;
	}

	public void setTimestamp(String timestamp) {
		this.timestamp = timestamp;
	}

	public String toHTML() {
		StringBuilder sb = new StringBuilder();

		sb.append(ResultConverter.getHTMLHeader());
		for (Result d : data) {
			sb.append(ResultConverter.toHTML(d));
		}

		sb.append(ResultConverter.getHTMLFooter());
		return sb.toString();
	}

	public String toXML() throws IOException {

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		PrintStream ps = new PrintStream(baos);

		try {
			Marshaller m = jaxbContext.createMarshaller();
			m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
			m.marshal(this, ps);
		} catch (JAXBException e) {
			throw new IOException(e);
		}

		return baos.toString();

	}

}
