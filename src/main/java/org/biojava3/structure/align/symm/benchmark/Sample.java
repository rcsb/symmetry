/*
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
 * Created on 2013-02-22
 *
 */
package org.biojava3.structure.align.symm.benchmark;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

/**
 * A collection of benchmark {@link Case Cases}.
 * @author dmyerstu
 */
@XmlRootElement(name = "CensusResults", namespace = "http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class Sample implements Serializable {

	private static final long serialVersionUID = 167984771894684447L;

	static JAXBContext jaxbContext;
	private List<Case> data;

	static {
		try {
			jaxbContext = JAXBContext.newInstance(Sample.class);
		} catch (Exception e) {
			throw new RuntimeException(e); // fatal
		}
	}

	public static Sample fromXML(File file) throws IOException {

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();
			FileInputStream fis = new FileInputStream(file);
			Sample sample = (Sample) un.unmarshal(fis);

			// due to a side effect by JAXB
			List<Case> newData = new ArrayList<Case>(sample.size());
			for (Case c : sample.getData()) {
				if (c != null) newData.add(c);
			}
			sample.setData(newData);

			return sample;

		} catch (JAXBException e) {
			throw new IOException(e);
		}

	}

	public static Sample fromXML(File[] files) throws IOException {
		Sample sample = new Sample();
		for (File file : files) {
			sample.getData().addAll(fromXML(file).getData());
		}
		return sample;
	}

	public Sample() {
		this.data = new ArrayList<Case>();
	}

	public Sample(List<Case> cases) {
		this.data = cases;
	}

	public boolean add(Case e) {
		return data.add(e);
	}

	public boolean addAll(Collection<? extends Case> c) {
		return data.addAll(c);
	}

	public boolean contains(Object o) {
		return data.contains(o);
	}

	public boolean containsAll(Collection<?> c) {
		return data.containsAll(c);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		Sample other = (Sample) obj;
		if (data == null) {
			if (other.data != null) return false;
		} else if (!data.equals(other.data)) return false;
		return true;
	}

	public List<Case> getData() {
		return data;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (data == null ? 0 : data.hashCode());
		return result;
	}

	public boolean remove(Object o) {
		return data.remove(o);
	}

	public boolean removeAll(Collection<?> c) {
		return data.removeAll(c);
	}

	public void setData(List<Case> cases) {
		this.data = cases;
	}

	public int size() {
		return data.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		int nAgree = 0;
		for (Case aCase : data) {
			final boolean agree = aCase.getOrder() == aCase.getKnownOrder();
			sb.append(aCase.getProtodomain() + ": " + (agree ? "agree" : "disagree")
					+ System.getProperty("line.seperator"));
			if (agree) nAgree++;
		}
		final double percentAgree = (double) nAgree / (double) data.size() * 100.0;
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(3);
		sb.append("Overall: " + nf.format(percentAgree) + " agree");
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
