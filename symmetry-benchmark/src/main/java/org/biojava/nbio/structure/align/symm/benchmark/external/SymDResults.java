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
 * Created on 2013-03-08
 *
 */
package org.biojava.nbio.structure.align.symm.benchmark.external;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

import org.biojava.nbio.structure.align.symm.census2.Alignment;
import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.symm.census2.Results;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Results of SymD.
 * 
 * @author dmyerstu
 * 
 */
@XmlRootElement(name = "CensusResults", namespace = "http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class SymDResults extends Results {

	private static JAXBContext jaxbContext;

	private static final Logger logger = LoggerFactory.getLogger( SymDResults.class );

	private static final long serialVersionUID = -6877332751979209323L;

	static {
		try {
			jaxbContext = JAXBContext.newInstance(SymDResults.class);
		} catch (Exception e) {
			throw new RuntimeException(e); // fatal
		}
	}

	/**
	 * Example output:
	 * 
	 * <pre>
	 * Program symd version 1.5b
	 * Number of residues read from the input file is 364.
	 * d1t3xa_  130 a.a. : Best(initial shift,N-aligned,T-score,Z-score)=( 106,   39,  18.998,   3.091)
	 * </pre>
	 * 
	 * @param output
	 * @return
	 */
	public static Result fromOutput13hw3(String output) {
		Result result = new Result();
		try {
			String[] lines = output.split("\n");
			String line = lines[lines.length - 1];
			String x = line.substring(line.lastIndexOf("(") + 1, line.length() - 1).trim();
			String[] values = x.split("[\\s,]+");
			result.setScopId(line.substring(0, line.indexOf(" ")));
			Alignment alignment = new Alignment();
			alignment.setInitialShift(Integer.parseInt(values[0]));
			alignment.setAlignLength(Integer.parseInt(values[1]));
			alignment.settScore(Float.parseFloat(values[2]));
			alignment.setSymDZScore(Float.parseFloat(values[3]));
			result.setAlignment(alignment);
		} catch (RuntimeException e) {
			throw new IllegalArgumentException("SymD returned strange output \"" + output + "\"", e);
		}
		return result;
	}

	/**
	 * Example output:
	 * 
	 * <pre>
	 * Program symd version 1.5b
	 * Number of residues read from the input file is 364.
	 * 1WOP  364 a.a. : Best(initial shift, N-aligned, N-non-self-aligned, Tm, Tmpr, Z1)=( 109,  140,  140,  134.07,  0.3683,  10.66)
	 * </pre>
	 * 
	 * @param output
	 * @return
	 */
	public static Result fromOutput15b(String output) {
		Result result = new Result();
		try {
			String[] lines = output.split("\n");
			String line = lines[lines.length - 1];
			String x = line.substring(line.lastIndexOf("(") + 1, line.length() - 1).trim();
			String[] values = x.split("[\\s,]+");
			result.setScopId(line.substring(0, line.indexOf(" ")));
			Alignment alignment = new Alignment();
			alignment.setInitialShift(Integer.parseInt(values[0]));
			alignment.setAlignLength(Integer.parseInt(values[1]));
			alignment.setnNonSelfAligned(Integer.parseInt(values[2]));
			alignment.settScore(Float.parseFloat(values[3]));
			alignment.setSymDTmScore(Float.parseFloat(values[4]));
			alignment.setSymDZScore(Float.parseFloat(values[5]));
			result.setAlignment(alignment);
		} catch (RuntimeException e) {
			throw new IllegalArgumentException("SymD returned strange output \"" + output + "\"", e);
		}
		return result;
	}

	public static SymDResults fromXML(File file) throws IOException {

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();
			FileInputStream fis = new FileInputStream(file);
			Results results = (Results) un.unmarshal(fis);

			// due to a side effect by JAXB
			List<Result> newData = new ArrayList<Result>(results.getData().size());
			for (Result result : results.getData()) {
				if (result != null)
					newData.add(result);
			}
			results.setData(newData);

			SymDResults symd = new SymDResults();
			for (Result result : results.getData())
				symd.add(result);

			return symd;

		} catch (JAXBException e) {
			throw new IOException(e);
		}

	}

	public static SymDResults fromXML(File[] files) throws IOException {
		SymDResults results = new SymDResults();
		for (File file : files) {
			results.getData().addAll(fromXML(file).getData());
		}
		return results;
	}

	public SymDResults() {
		super();
	}

	@Override
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
