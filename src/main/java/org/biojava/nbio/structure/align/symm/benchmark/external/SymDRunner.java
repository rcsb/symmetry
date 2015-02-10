package org.biojava.nbio.structure.align.symm.benchmark.external;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.benchmark.SampleBuilder;
import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.utils.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Builds {@link SymDResults} by running SymD on a list of domains.
 * 
 * @author dmyerstu
 */
public class SymDRunner {

	private static final Logger logger = LoggerFactory.getLogger( SymDRunner.class );

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

	public static float getTmScoreFromFastaFile(String scopId, Structure structure) throws SymDException {
		File fastaFile = new File(scopId + "-best.fasta");
		try {
			AFPChain afpChain = SymDFasta.getAlignment(fastaFile);
			if (afpChain == null)
				throw new SymDException("AFPChain is null");
			return (float) afpChain.getTMScore();
		} catch (Exception e) {
			throw new SymDException("The FASTA file was wrong or could not be converted to an AFPChain", e);
		}
	}

	public static void main(String[] args) {
		if (args.length > 3 || args.length > 4) {
			System.err.println("Usage: " + SymDResults.class.getSimpleName()
					+ " SymD-path names-file output-file [1.5b]");
			return;
		}
		final String symDPath = args[0];
		final File namesFile = new File(args[1]);
		final File outputFile = new File(args[2]);
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75);
		ScopFactory.setScopDatabase(scop);
		AtomCache cache = new AtomCache();
		boolean updated = false;
		if (args.length > 3) {
			if (args[3].equalsIgnoreCase("1.5b") || args[3].equalsIgnoreCase("updated")
					|| args[3].equalsIgnoreCase("true"))
				updated = true;
		}
		SymDRunner runner = new SymDRunner(cache, scop, symDPath, outputFile, updated);
		runner.run(namesFile);
	}

	/**
	 * Useful because some processes won't produce input streams that can be
	 * buffered (instead, they discard their output immediately).
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
				logger.warn(e.getMessage(),e);
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

	private int printFrequency = 20;
	private boolean isUpdated;
	private String symdPdbFilesDir;
	private AtomCache cache;

	private String symdPath;

	private File outputFile;

	private ScopDatabase scop;

	private double timeTaken = 0;

	public SymDRunner(AtomCache cache, ScopDatabase scop, String symdPath, File outputFile, boolean isUpdated) {
		super();
		this.isUpdated = isUpdated;
		this.cache = cache;
		this.symdPath = symdPath;
		symdPdbFilesDir = symdPdbFilesDir.substring(0, symdPdbFilesDir.lastIndexOf('/'));
		this.outputFile = outputFile;
		this.scop = scop;
	}

	public double getTimeTaken() {
		return timeTaken;
	}

	public void printResults(SymDResults results) {
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

	public void run(File namesFile) {

		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		List<String> names;
		try {
			names = SampleBuilder.getNames(namesFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		for (String d : names) {
			domains.add(scop.getDomainByScopID(d));
		}

		runOnDomains(domains);

	}

	public SymDResults runOnDomains(List<ScopDomain> scopDomains) {
		SymDResults results = getResultsFromPrevRun(outputFile);
		long timeTaken = 0;
		int i = 0;
		int nSuccess = 0;
		for (ScopDomain domain : scopDomains) {
			try {
				Structure structure;
				try {
					structure = cache.getStructure(domain.getScopId());
				} catch (StructureException e) {
					throw new RuntimeException("Could not get Structure for domain " + domain.getScopId(), e);
				} catch (IOException e) {
					throw new RuntimeException("Could not create PDB file for domain " + domain.getScopId(), e);
				}
				final File file = new File(symdPdbFilesDir + domain.getScopId() + ".pdb");
				if (!file.exists()) {
					try {
						BufferedWriter bw = new BufferedWriter(new FileWriter(file));
						bw.write(structure.toPDB());
						bw.close();
					} catch (IOException e) {
						throw new RuntimeException("Could not create PDB file for domain " + domain.getScopId(), e);
					}
				}
				Result result;
				try {
					long startTime = System.currentTimeMillis();
					if (isUpdated) {
						result = runSymD15b(domain, structure, file.getPath());
					} else {
						result = runSymD13hw3(domain, structure, file.getPath());
					}
					long endTime = System.currentTimeMillis();
					timeTaken += endTime - startTime;
					nSuccess++;
				} catch (SymDException e) {
					logger.error("SymD failed on " + domain.getScopId(), e);
					continue;
				}
				results.add(result);
			} finally {
				i++;
			}
			if (i % printFrequency == 0) {
				printResults(results);
			}
		}
		printResults(results);
		this.timeTaken = (double) timeTaken / (double) nSuccess;
		results.setMeanSecondsTaken(this.timeTaken);
		return results;
	}

	public Result runSymD13hw3(ScopDomain domain, Structure structure, String pdbFilePath) throws SymDException {
		final String[] cmd = new String[] { symdPath, pdbFilePath };
		final String output = runCmd(cmd); // waits for completion
		try {
			return SymDResults.fromOutput13hw3(output);
		} catch (IllegalArgumentException e) {
			throw new SymDException("SymD failed on " + pdbFilePath, e);
		}
	}

	private Result runSymD15b(ScopDomain domain, Structure structure, String pdbFilePath) throws SymDException {
		final String[] cmd = new String[] { symdPath, pdbFilePath };
		final String output = runCmd(cmd); // waits for completion
		Result result;
		try {
			result = SymDResults.fromOutput15b(output);
		} catch (IllegalArgumentException e) {
			throw new SymDException("SymD failed on " + pdbFilePath, e);
		}
		try {
			float tmScore = getTmScoreFromFastaFile(domain.getScopId(), structure);
			result.getAlignment().setTmScore(tmScore);
		} catch (SymDException e) {
			logger.error("Couldn't set TM-score for " + domain.getScopId(), e);
		}
		return result;
	}

	public void setPrintFrequency(int printFrequency) {
		this.printFrequency = printFrequency;
	}

}
