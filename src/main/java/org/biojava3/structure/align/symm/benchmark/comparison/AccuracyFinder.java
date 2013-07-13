/*
 * BioJava development code
 * 
 * This code may be freely distributed and modified under the terms of the GNU Lesser General Public Licence. This
 * should be distributed with the code. If you do not have a copy, see:
 * 
 * http://www.gnu.org/copyleft/lesser.html
 * 
 * Copyright for this code is held jointly by the individual authors. These should be listed in @author doc comments.
 * 
 * For more information on the BioJava project and its aims, or to join the biojava-l mailing list, visit the home page
 * at:
 * 
 * http://www.biojava.org/
 * 
 * Created on 2013-02-22
 */
package org.biojava3.structure.align.symm.benchmark.comparison;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava3.structure.align.symm.benchmark.Case;
import org.biojava3.structure.align.symm.benchmark.Sample;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

/**
 * Determines the false-positive, true-positive, false-negative, and true-negative frequencies of a benchmark
 * {@link Sample} using a particular Significance object (see symmetry project).
 * Considers only the binary choice of symmetric versus asymmetric.
 * @author dmyerstu
 * @see OrderAccuracy, which determines the accuracy of determining order rather than just symmetric/asymmetric
 */
public class AccuracyFinder {

	private static final Logger logger = LogManager.getLogger(AccuracyFinder.class.getName());

	public static void main(String[] args) throws IOException {
		if (args.length != 1 && args.length != 2) {
			System.err.println("Usage: " + AccuracyFinder.class.getSimpleName() + " input-file [significance-class] [significance-method]");
			return;
		}
		File input = new File(args[0]);
		Significance sig = SignificanceFactory.rotationallySymmetricSmart();
		if (args.length > 2) {
			logger.info("Using alternate Significance "  + args[1] + "." + args[2]);
			sig = SignificanceFactory.fromMethod(args[1], args[2]);
		} else if (args.length > 1) {

			logger.info("Using alternate Significance "  + args[1]);
			sig = SignificanceFactory.fromMethod(SignificanceFactory.class.getName(), args[1]);
		}
		AccuracyFinder finder = new AccuracyFinder(input, sig);
		System.out.println(finder);
	}

	public AccuracyFinder(File input, Significance sig) throws IOException {
		this(Sample.fromXML(input), sig);
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

	public AccuracyFinder(Sample sample, Significance sig) {
		for (Case c : sample.getData()) {
			try {
				if (c.getKnownInfo().hasRotationalSymmetry()) {
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
		String tps = nf.format((double) tp / (double) (tp + tn + fp + fn) * 100.0) + "%";
		String fps = nf.format((double) fp / (double) (tp + tn + fp + fn) * 100.0) + "%";
		String tns = nf.format((double) tn / (double) (tp + tn + fp + fn) * 100.0) + "%";
		String fns = nf.format((double) fn / (double) (tp + tn + fp + fn) * 100.0) + "%";
		sb.append("True positives: \t" + tp + "\t(" + tps + ")\n");
		sb.append("True negatives: \t" + tn + "\t(" + tns + ")\n");
		sb.append("False positives: \t" + fp + "\t(" + fps + ")\n");
		sb.append("False negatives: \t" + fn + "\t(" + fns + ")\n");
		return sb.toString();
	}

}
