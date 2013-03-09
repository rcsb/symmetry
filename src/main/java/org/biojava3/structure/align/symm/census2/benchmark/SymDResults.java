package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
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

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Alignment;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;

@XmlRootElement(name = "SymDResults", namespace = "http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class SymDResults extends Results {

	private static JAXBContext jaxbContext;

	private static final long serialVersionUID = -6877332751979209323L;
	static final Logger logger = Logger.getLogger(SymDResults.class.getPackage().getName());

	static {
		try {
			jaxbContext = JAXBContext.newInstance(SymDResults.class);
		} catch (Exception e) {
			throw new RuntimeException(e); // fatal
		}
	}

	static {
		BasicConfigurator.configure();
		logger.setLevel(Level.DEBUG);
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
			String line = output.split("\n")[2];
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
		Census.setBerkeleyScop(pdbDir);
		AtomCache cache = new AtomCache(pdbDir, false);
		writeToFile(symDPath, outputFile, cache, namesFile);
	}

	public static SymDResults runSymD(String symDPath, String pdbFilesPath, AtomCache cache, List<ScopDomain> scopDomains) {
		SymDResults results = new SymDResults();
		if (!pdbFilesPath.endsWith("/"))
			pdbFilesPath += "/";
		for (ScopDomain domain : scopDomains) {
			final File file = new File(pdbFilesPath + domain.getScopId() + ".pdb");
			//if (!file.exists()) { TODO why is this breaking everything?
				try {
					Structure structure = cache.getStructure(domain.getScopId());
					BufferedWriter bw = new BufferedWriter(new FileWriter(file));
					bw.write(structure.toPDB());
					bw.close();
				} catch (IOException e) {
					throw new RuntimeException("Could not create PDB file for domain " + domain.getScopId(), e);
				} catch (StructureException e) {
					throw new RuntimeException("Could not get Structure for domain " + domain.getScopId(), e);
				}
			//}
			Result result;
			try {
				result = runSymD(symDPath, file.getPath());
			} catch (SymDException e) {
				logger.error("SymD failed on " + domain.getScopId(), e);
				continue;
			}
			results.add(result);
		}
		return results;
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
		try {
			BufferedReader br = new BufferedReader(new FileReader(lineByLine));
			String line = "";
			while ((line = br.readLine()) != null) {
				if (line.trim().isEmpty())
					continue;
				ScopDomain domain = scop.getDomainByScopID(line);
				if (domain == null) {
					logger.error("No SCOP domain with id " + line + " was found");
				} else {
					domains.add(domain);
				}
			}
			br.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		writeToFile(symDPath, domains, cache, outputFile);

	}

	public static void writeToFile(String symDPath, List<ScopDomain> scopDomains, AtomCache cache, File outputFile) {
		final String filesPath = symDPath.substring(0, symDPath.lastIndexOf('/'));
		SymDResults results = SymDResults.runSymD(symDPath, filesPath, cache, scopDomains);
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
		Thread drainer = new Thread() {
			@Override
			public void run() {
				int c;
				do {
					try {
						c = out.read();
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
					if (c >= 0)
						sb.append((char) c);
				} while (c >= 0);
			}
		};
		drainer.start();
		while (true) {
			try {
				process.waitFor();
				break;
			} catch (InterruptedException e) {
			}
		}
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
