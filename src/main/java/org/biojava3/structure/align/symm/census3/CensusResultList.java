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
package org.biojava3.structure.align.symm.census3;

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
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * A census of symmetry; that is, a list of {@link CensusResult CensusResults}.
 * @author dmyersturnbull
 */
@XmlRootElement(name = "CensusResultList", namespace = "http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class CensusResultList implements Serializable {

	private static final Logger logger = LogManager.getLogger(CensusResultList.class.getSimpleName());

	private static final long serialVersionUID = -5517546595033480440L;
	private static JAXBContext jaxbContext;
	
	private List<CensusResult> entries = new ArrayList<CensusResult>();

	private String startingTime;
	private float meanSecondsTaken;

	static {
		try {
			jaxbContext = JAXBContext.newInstance(CensusResultList.class);
		} catch (Exception e) {
			throw new RuntimeException(e); // fatal
		}
	}

	public static CensusResultList fromXML(File file) throws IOException {

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();
			FileInputStream fis = new FileInputStream(file);
			CensusResultList results = (CensusResultList) un.unmarshal(fis);

			// due to a side effect by JAXB
			List<CensusResult> newData = new ArrayList<CensusResult>(results.getEntries().size());
			
			for (CensusResult result : results.getEntries()) {
				if (result != null) newData.add(result);
			}
			results.setEntries(newData);

			return results;

		} catch (JAXBException e) {
			throw new IOException(e);
		}

	}
	
	public static CensusResultList fromXML(File[] files) throws IOException {
		CensusResultList results = new CensusResultList();
		// don't keep duplicate SCOP Ids
		Set<String> scopIds = new HashSet<String>();
		float meanSecondsTaken = 0f;
		for (File file : files) {
			CensusResultList old = fromXML(file);
			logger.debug("Taking " + old.size() + " results from " + file.getName());
			meanSecondsTaken += old.getMeanSecondsTaken();
			for (CensusResult result : old.getEntries()) {
				if (!scopIds.contains(result.getId())) {
					scopIds.add(result.getId());
					results.add(result);
				} else {
					logger.warn("Found duplicate " + result.getId() + " in file " + file.getName());
				}
			}
		}
		results.setMeanSecondsTaken(meanSecondsTaken / files.length);
		return results;
	}

	public static CensusResultList getExistingResults(File file) throws IOException {
		if (file.exists()) {
			return CensusResultList.fromXML(file);
		}
		return null;
	}

	public static CensusResultList getExistingResults(String file) throws IOException {
		return getExistingResults(new File(file));
	}

	public CensusResultList() {
		startingTime = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date());
	}

	public boolean add(CensusResult e) {
		return entries.add(e);
	}

	public boolean addAll(Collection<? extends CensusResult> c) {
		return entries.addAll(c);
	}

	public int size() {
		return entries.size();
	}

	public boolean contains(Object o) {
		return entries.contains(o);
	}

	public boolean containsAll(Collection<?> c) {
		return entries.containsAll(c);
	}

	@XmlElement(name = "entry")
	public List<CensusResult> getEntries() {
		return entries;
	}

	public String getStartingTime() {
		return startingTime;
	}

	public void setEntries(List<CensusResult> entries) {
		this.entries = entries;
	}

	public void setStartingTime(String timestamp) {
		this.startingTime = timestamp;
	}

	public float getMeanSecondsTaken() {
		return meanSecondsTaken;
	}

	public void setMeanSecondsTaken(float meanSecondsTaken) {
		this.meanSecondsTaken = meanSecondsTaken;
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
