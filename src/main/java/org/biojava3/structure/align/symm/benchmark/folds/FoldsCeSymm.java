package org.biojava3.structure.align.symm.benchmark.folds;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopInstallation;
import org.biojava3.structure.align.symm.census2.NamesCensus;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

/**
 * Runs CE-Symm on a list of domains specified in the SymD paper (Kim et. al.).
 * The files are in {@code src/main/resources/Guerler_folds}.
 * @author dmyerstu
 */
public class FoldsCeSymm {

	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + FoldsCeSymm.class.getSimpleName() + " output-dir");
			return;
		}
		run(args[0]);
	}
	
	private static final String RESOURCE_DIR = "src/main/resources/Guerler_folds/";
	
	public static void run(String dir) throws IOException {
		ScopInstallation scop = new ScopInstallation();
		scop.setScopVersion("1.73");
		ScopFactory.setScopDatabase(scop);
		Significance tm = SignificanceFactory.forCeSymmTm();
		Significance ord = SignificanceFactory.forCeSymmOrd();
		new File(dir).mkdirs();
		String[] folds = new String[] {"a.24", "b.1", "b.11", "b.42", "b.69", "c.1", "d.131", "d.58"};
		System.out.println("fold\tTM\tord\tN");
		for (String fold : folds) {
			File censusFile = new File(dir + fold + ".xml");
			File lineByLine = new File(RESOURCE_DIR + fold + "_names.list");
			NamesCensus.buildDefault(censusFile, lineByLine, false);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir + "stats.txt"), true)));
			int x = 0;
			int y = 0;
			Results results = Results.fromXML(censusFile);
			for (Result result : results.getData()) {
				if (tm.isSignificant(result)) x++;
				if (ord.isSignificant(result)) y++;
			}
			pw.println(fold + "\t" + x + "\t" + y + "\t" + results.size());
			pw.close();
		}
	}
	
}
