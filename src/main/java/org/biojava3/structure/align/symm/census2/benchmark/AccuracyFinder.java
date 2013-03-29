package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

public class AccuracyFinder {

	static final Logger logger = Logger.getLogger(Case.class.getPackage().getName());

	static {
		BasicConfigurator.configure();
		logger.setLevel(Level.DEBUG);
	}

	public static void main(String[] args) throws IOException {
		AccuracyFinder finder = new AccuracyFinder(new File(args[0]));
		System.out.println(finder);
	}

	public AccuracyFinder(File input) throws IOException {
		this(Sample.fromXML(input));
	}

	public int getTp() {
		return tp;
	}

	public int getFn() {
		return fn;
	}

	public int getFp() {
		return fp;
	}

	public int getTn() {
		return tn;
	}

	private int tp = 0;
	private int fn = 0;
	private int fp = 0;
	private int tn = 0;
	
	private static Significance sig;
	static {
		sig = new Significance() {
			private static final double cutoff = 0.55;
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return true; // whatever
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				return true; // whatever
			}

			@Override
			public boolean isSignificant(Result result) {
				if (result.getOrder() == null || result.getAlignment() == null || result.getAxis() == null) return false;
				return result.getAlignment().getTmScore() >= cutoff && result.getOrder() > 1;
			}
		};
	}
	
	
	public AccuracyFinder(Sample sample) {
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
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		NumberFormat nf = new DecimalFormat();
		nf.setMaximumFractionDigits(2);
		String tps = nf.format((double) tp / (double) (tp+tn+fp+fn) * 100.0) + "%";
		String fps = nf.format((double) fp / (double) (tp+tn+fp+fn) * 100.0) + "%";
		String tns = nf.format((double) tn / (double) (tp+tn+fp+fn) * 100.0) + "%";
		String fns = nf.format((double) fn / (double) (tp+tn+fp+fn) * 100.0) + "%";
		sb.append("True positives: \t" + tp + "\t(" + tps + ")\n");
		sb.append("True negatives: \t" + tn + "\t(" + tns + ")\n");
		sb.append("False positives: \t" + fp + "\t(" + fps + ")\n");
		sb.append("False negatives: \t" + fn + "\t(" + fns + ")\n");
		return sb.toString();
	}

}
