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
package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FastaAFPChainConverter;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Alignment;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.utils.FileUtils;

/**
 * Results of SymD.
 * @author dmyerstu
 *
 */
@XmlRootElement(name = "CensusResults", namespace = "http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class SymDResults extends Results {

	private static JAXBContext jaxbContext;

	private static final Logger logger = LogManager.getLogger(SymDResults.class.getName());

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
	 * 1WOP  364 a.a. : Best(initial shift, N-aligned, N-non-self-aligned, Tm, Tmpr, Z1)=( 109,  140,  140,  134.07,  0.3683,  10.66)
	 * </pre>
	 * 
	 * @param output
	 * @return
	 */
	public static Result fromOutput(String output) {
		Result result = new Result();
		try {
			String[] lines = output.split("\n");
			String line = lines[lines.length-1];
			String x = line.substring(line.lastIndexOf("(") + 1, line.length() - 1).trim();
			String[] values = x.split("[\\s,]+");
			result.setScopId(line.substring(0, line.indexOf(" ")));
			Alignment alignment = new Alignment();
			alignment.setInitialShift(Integer.parseInt(values[0]));
			alignment.setAlignLength(Integer.parseInt(values[1]));
			alignment.setnNonSelfAligned(Integer.parseInt(values[2]));
			alignment.setAlternateTm(Float.parseFloat(values[3]));
			alignment.setTmpr(Float.parseFloat(values[4]));
			alignment.setzScore(Float.parseFloat(values[5]));
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
			SymDResults results = (SymDResults) un.unmarshal(fis);

			// due to a side effect by JAXB
			List<Result> newData = new ArrayList<Result>(results.getData().size());
			for (Result result : results.getData()) {
				if (result != null)
					newData.add(result);
			}
			results.setData(newData);

			return results;

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

	public static void main(String[] args) {
		final String pdbDir = args[0];
		final String symDPath = args[1];
		final File namesFile = new File(args[2]);
		final File outputFile = new File(args[3]);
		ScopFactory.setScopDatabase(new BerkeleyScopInstallation());
		AtomCache cache = new AtomCache(pdbDir, false);
		writeToFile(symDPath, namesFile, cache, outputFile);
	}

	private static SymDResults getResultsFromPrevRun(File file) {
		if (file.exists() && file.length() > 0) {
			try {
				SymDResults results = SymDResults.fromXML(file);
				logger.info("Found " + results.size() + " previous results from " + file.getPath());
				return results;
			} catch (IOException e) {
				final Date date = new Date();
				try {
					logger.warn("Could not load file " + file.getPath() + ". Starting from scratch.", e);
					FileUtils.copy(file, new File(file.getPath() + " __backup " + date));
					file.delete();
				} catch (IOException e1) {
					throw new RuntimeException("Could not read census file, and could not backup previous file", e1);
				}
			}
		}
		return new SymDResults();
	}

	public static SymDResults runSymD(String symDPath, String pdbFilesPath, AtomCache cache, List<ScopDomain> scopDomains, File resultsFile) {
		SymDResults results = getResultsFromPrevRun(resultsFile);
		long timeTaken = 0;
		int nSuccess = 0;
		if (!pdbFilesPath.endsWith("/"))
			pdbFilesPath += "/";
		for (ScopDomain domain : scopDomains) {
			Structure structure = null;
			final File file = new File(pdbFilesPath + domain.getScopId() + ".pdb");
			if (!file.exists()) {
				try {
					structure = cache.getStructure(domain.getScopId());
					BufferedWriter bw = new BufferedWriter(new FileWriter(file));
					bw.write(structure.toPDB());
					bw.close();
				} catch (IOException e) {
					throw new RuntimeException("Could not create PDB file for domain " + domain.getScopId(), e);
				} catch (StructureException e) {
					throw new RuntimeException("Could not get Structure for domain " + domain.getScopId(), e);
				}
			}
			Result result;
			try {
				long startTime = System.currentTimeMillis();
				result = runSymD(symDPath, file.getPath());
				float tmScore = getTmScoreFromFastaFile(symDPath, domain.getScopId(), structure);
				result.getAlignment().setTmScore(tmScore);
				long endTime = System.currentTimeMillis();
				timeTaken += (endTime - startTime);
				nSuccess++;
			} catch (SymDException e) {
				logger.error("SymD failed on " + domain.getScopId(), e);
				continue;
			}
			results.add(result);
		}
		System.err.println("AVG TIME TAKEN: " + ((double) timeTaken / (double) nSuccess));
		return results;
	}

	public static float getTmScoreFromFastaFile(String symDPath, String scopId, Structure structure) throws SymDException {
		File fastaFile = new File(new File(symDPath).getParent() + scopId + "-best.fasta");
		try {
			AFPChain afpChain = FastaAFPChainConverter.fastaFileToAfpChain(fastaFile, structure, structure);
			return (float) afpChain.getTMScore();
		} catch (Exception e) {
			throw new SymDException("The FASTA file was wrong or could not be converted to an AFPChain");
		}
	}
	
	public static Result runSymD(String symDPath, String pdbFilePath) throws SymDException {
		final String[] cmd = new String[] {symDPath, pdbFilePath};
		final String output = runCmd(cmd); // waits for completion
		try {
			return fromOutput(output);
		} catch (IllegalArgumentException e) {
			throw new SymDException("SymD failed on " + pdbFilePath, e);
		}
	}

	public static void writeToFile(String symDPath, File lineByLine, AtomCache cache, File outputFile) {

		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		ScopDatabase scop = ScopFactory.getSCOP();
		Map<String,KnownInfo> infos;
		try {
			infos = SampleBuilder.getOrders(lineByLine);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		for (String d : infos.keySet()) {
			domains.add(scop.getDomainByScopID(d));
		}

		writeToFile(symDPath, domains, cache, outputFile);

	}

	public static void writeToFile(String symDPath, List<ScopDomain> scopDomains, AtomCache cache, File outputFile) {
		final String filesPath = symDPath.substring(0, symDPath.lastIndexOf('/'));
		SymDResults results = SymDResults.runSymD(symDPath, filesPath, cache, scopDomains, outputFile);
		String xml;
		try {
			xml = results.toXML();
		} catch (IOException e) {
			throw new RuntimeException("Couldn't get XML results", e);
		}
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(outputFile));
			bw.write(xml);
			bw.close();
		} catch (IOException e) {
			throw new RuntimeException("Couldn't write XML results to " + outputFile, e);
		}
	}

	/**
	 * Useful because some processes won't produce input streams that can be buffered (instead, they discard their output immediately).
	 * 
	 * @param cmd
	 * @return
	 */
	private static String runCmd(String[] cmd) {
		Process process;
		try {
			process = new ProcessBuilder(cmd).redirectErrorStream(true).start();
		} catch (IOException e) {
			throw new RuntimeException("Could not create execution process", e);
		}
		final InputStream out = process.getInputStream();
		final StringBuilder sb = new StringBuilder();
		while (true) {
			try {
				process.waitFor();
				break;
			} catch (InterruptedException e) {
			}
		}

		int c;
		do {
			try {
				c = out.read();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			if (c >= 0)
				sb.append((char) c);
		} while (c != -1);
		return sb.toString();
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
