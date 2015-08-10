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
package org.biojava.nbio.structure.align.symm.census2;

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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

import org.biojava.nbio.structure.align.symm.census2.utils.ResultConverter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A census (list of census {@link Result Results}).
 * @author dmyersturnbull
 * @deprecated
 */
@Deprecated
@XmlRootElement(name = "CensusResults", namespace = "http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class Results implements Serializable {

	private final static Logger logger = LoggerFactory.getLogger(Results.class);

	private static final long serialVersionUID = -5517546595033480440L;
	private static JAXBContext jaxbContext;
	private List<Result> data = new ArrayList<Result>();

	private String timestamp;
	private double meanSecondsTaken;

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
			results.add(new SimpleResult(r));
		}
		return results;
	}

	public static Results fromXML(File[] files) throws IOException {
		Results results = new Results();
		// don't keep duplicate SCOP Ids
		Set<String> scopIds = new HashSet<String>();
		double meanSecondsTaken = 0;
		for (File file : files) {
			Results old = fromXML(file);
			logger.debug("Taking " + old.size() + " results from " + file.getName());
			meanSecondsTaken += old.getMeanSecondsTaken();
			for (Result result : old.getData()) {
				if (!scopIds.contains(result.getScopId())) {
					scopIds.add(result.getScopId());
					results.add(result);
				} else {
					logger.warn("Found duplicate " + result.getScopId() + " in file " + file.getName());
				}
			}
		}
		results.setMeanSecondsTaken(meanSecondsTaken / files.length);
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

	public double getMeanSecondsTaken() {
		return meanSecondsTaken;
	}

	public void setMeanSecondsTaken(double meanSecondsTaken) {
		this.meanSecondsTaken = meanSecondsTaken;
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

	public static Results fromXML(String queryFile) throws IOException {
		return fromXML(new File(queryFile));
	}

}
