package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

public class AccuracyFinder {

	static final Logger logger = Logger.getLogger(Case.class.getPackage().getName());

	static {
		BasicConfigurator.configure();
		logger.setLevel(Level.DEBUG);
	}

	public static void main(String[] args) throws IOException {
		findAccuracy(new File(args[0]), System.out);
	}

	public static void findAccuracy(File input, PrintStream ps) throws IOException {
		findAccuracy(Sample.fromXML(input), ps);
	}

	public static void findAccuracy(Sample sample, PrintStream ps) {
		int tp = 0, fp = 0, tn = 0, fn = 0;
		final Significance sig = SignificanceFactory.getGenerallySymmetric();
		for (Case c : sample.getData()) {
			try {
				if (c.hasKnownSymmetry()) {
					if (sig.isSignificant(c.getResult())) {
						tp++;
					} else {
						fn++;
					}
				} else {
					if (sig.isSignificant(c.getResult())) {
						fp++;
					} else {
						tn++;
					}
				}
			} catch (RuntimeException e) {
				logger.error("Encountered an error on " + c.getScopId(), e);
			}
		}
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(2);
		String tps = nf.format((double) tp / (double) (tp+tn+fp+fn) * 100.0) + "%";
		String fps = nf.format((double) fp / (double) (tp+tn+fp+fn) * 100.0) + "%";
		String tns = nf.format((double) tn / (double) (tp+tn+fp+fn) * 100.0) + "%";
		String fns = nf.format((double) fn / (double) (tp+tn+fp+fn) * 100.0) + "%";
		ps.println("True positives: \t" + tp + "\t(" + tps + ")");
		ps.println("True negatives: \t" + tn + "\t(" + tns + ")");
		ps.println("False positives: \t" + fp + "\t(" + fps + ")");
		ps.println("False negatives: \t" + fn + "\t(" + fns + ")");
	}

}
