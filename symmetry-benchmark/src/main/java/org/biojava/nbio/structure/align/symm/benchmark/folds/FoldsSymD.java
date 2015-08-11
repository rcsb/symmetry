package org.biojava.nbio.structure.align.symm.benchmark.folds;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.nbio.structure.align.symm.benchmark.external.SymDRunner;
import org.biojava.nbio.structure.align.symm.census2.NamesCensus;
import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.symm.census2.Results;
import org.biojava.nbio.structure.align.symm.census2.Significance;
import org.biojava.nbio.structure.align.symm.census2.SignificanceFactory;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.scop.ScopInstallation;

/**
 * Runs SymD on a list of domains specified in the SymD paper (Kim et. al.).
 * The files are in {@code src/main/resources/Guerler_folds}.
 * @author dmyerstu
 */
public class FoldsSymD {

	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: " + FoldsSymD.class.getSimpleName() + " output-dir symd-path");
			return;
		}
		run(args[0], args[1]);
	}

	private static final String RESOURCE_DIR = "src/main/resources/Guerler_folds/";
	
	public static void run(String dir, String symdPath) throws IOException {
		final AtomCache cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);
		ScopInstallation scop = new ScopInstallation();
		scop.setScopVersion("1.73");
		ScopFactory.setScopDatabase(scop);
		Significance sig = SignificanceFactory.generallySymmetric();
		new File(dir).mkdirs();
		String[] folds = new String[] {"a.24", "b.1", "b.11", "b.42", "b.69", "c.1", "d.131", "d.58"};
		for (String fold : folds) {
			File censusFile = new File(dir + fold + ".xml");
			File lineByLine = new File(RESOURCE_DIR + fold + "_names.list");
			SymDRunner runner = new SymDRunner(cache, scop, symdPath, censusFile, true);
			runner.run(lineByLine);
			NamesCensus.buildDefault(censusFile, lineByLine, false);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir + "stats.txt"), true)));
			int x = 0;
			Results results = Results.fromXML(censusFile);
			for (Result result : results.getData()) {
				if (sig.isSignificant(result)) x++;
			}
			pw.println(fold + "\t" + x + "\t" + results.size());
			pw.close();
		}
	}

}
