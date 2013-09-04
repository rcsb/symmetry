package demo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopInstallation;
import org.biojava3.structure.align.symm.benchmark.external.SymDRunner;
import org.biojava3.structure.align.symm.census2.NamesCensus;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

/**
 * Runs SymD on a list of domains specified in the SymD paper (Kim et. al.).
 * The files are in {@code src/main/resources/Guerler_folds}.
 * @author dmyerstu
 */
public class GuerlerFoldsSymD {

	private static String DIR = "/Users/dmyerstu/Desktop/folds/";
	private static String RESOURCE_DIR = "/Users/dmyerstu/Documents/workspace/symmetry-benchmark/src/main/resources/Guerler_folds/";
	private static String SYMD_RESOURCE_PATH = "/Users/dmyerstu/Desktop/SymD/symd_result_files/";
	private static String SYMD_PATH = "/Users/dmyerstu/Desktop/SymD/symd_revised";
	private static File STATS_FILE = new File("/Users/dmyerstu/Desktop/folds/stats.txt");
	
	public static void main(String[] args) throws IOException {
		final AtomCache cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);
		ScopInstallation scop = new ScopInstallation();
		scop.setScopVersion("1.73");
		ScopFactory.setScopDatabase(scop);
		Significance sig = SignificanceFactory.generallySymmetric();
		new File(DIR).mkdirs();
		String[] folds = new String[] {"a.24", "b.1", "b.11", "b.42", "b.69", "c.1", "d.131", "d.58"};
		for (String fold : folds) {
			File censusFile = new File(DIR + fold + ".xml");
			File lineByLine = new File(RESOURCE_DIR + fold + "_names.list");
			SymDRunner runner = new SymDRunner(cache, scop, SYMD_PATH, censusFile, true);
			runner.run(lineByLine);
			NamesCensus.buildDefault(censusFile, lineByLine, false);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(STATS_FILE, true)));
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
