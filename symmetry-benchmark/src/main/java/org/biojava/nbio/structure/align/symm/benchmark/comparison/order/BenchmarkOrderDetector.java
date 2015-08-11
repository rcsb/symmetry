package org.biojava.nbio.structure.align.symm.benchmark.comparison.order;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.benchmark.Case;
import org.biojava.nbio.structure.align.symm.benchmark.Sample;
import org.biojava.nbio.structure.align.symm.census3.CensusResult;
import org.biojava.nbio.structure.align.symm.census3.CensusSignificanceFactory;
import org.biojava.nbio.structure.align.symm.census3.CensusSymmetryGroup;
import org.biojava.nbio.structure.align.symm.order.AngleOrderDetectorPlus;
import org.biojava.nbio.structure.align.symm.order.HybridOrderDetector;
import org.biojava.nbio.structure.align.symm.order.MultipassOrderDetector;
import org.biojava.nbio.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava.nbio.structure.align.symm.order.OrderDetector;
import org.biojava.nbio.structure.align.symm.order.PeakCountingOrderDetector;
import org.biojava.nbio.structure.align.symm.order.RotationOrderDetector;
import org.biojava.nbio.structure.align.symm.order.RotationOrderDetector.RotationOrderMethod;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Test the performance of a {@link OrderDetector} against a benchmark.
 * @author dmyersturnbull
 */
public class BenchmarkOrderDetector {
	private static final Logger logger = LoggerFactory.getLogger(BenchmarkOrderDetector.class);

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length >2) {
			System.err.println("Usage: " + BenchmarkOrderDetector.class.getSimpleName() + " input-benchmark-file [output-file]");
			return;
		}
		String inputFilename = args[0];
		String outputFilename = null;
		if(args.length>1) {
			outputFilename = args[1];
		}

		final int maxOrder =8;

		List<OrderDetector> detectors = new LinkedList<OrderDetector>();
		for(RotationOrderMethod method : RotationOrderMethod.values()) {
			detectors.add( new RotationOrderDetector(maxOrder,method) );
		}
//		detectors.add( new AngleOrderDetector(maxOrder,2*Math.PI/16) );
		detectors.add( new AngleOrderDetectorPlus(maxOrder,2*Math.PI/16,false) ); //should be equivalent
//		detectors.add( new AngleOrderDetectorPlus(maxOrder,2*Math.PI,false) ); // infinite error
//		detectors.add( new AngleOrderDetector(maxOrder,Math.toRadians(5)) );
		detectors.add( new AngleOrderDetectorPlus(maxOrder,Math.toRadians(5),false) );
//		detectors.add( new AngleOrderDetector(maxOrder,Math.toRadians(1)) );
		detectors.add( new AngleOrderDetectorPlus(maxOrder,Math.toRadians(1),false) );
//		detectors.add( new AngleOrderDetectorPlus(maxOrder,1,true));// infinite error
		detectors.add( new AngleOrderDetectorPlus(maxOrder,1/9.,true)); // should be equivalent
		detectors.add( new AngleOrderDetectorPlus(maxOrder,1/15.,true));
		detectors.add( new SequenceFunctionOrderDetector() );
		detectors.add( new PeakCountingOrderDetector(8) );
		detectors.add( new HybridOrderDetector(maxOrder, Math.PI/16, false, .85));
		detectors.add( new HybridOrderDetector(maxOrder, .07, false, .85));
		detectors.add( new HybridOrderDetector(maxOrder, .07, false, .6));
		detectors.add( new HybridOrderDetector(maxOrder, .07, false, 1.));
		detectors.add( new HybridOrderDetector(maxOrder, .05, false, .85));
		
		detectors.add(new MultipassOrderDetector(maxOrder));

		// Read in benchmark
		AtomCache cache = new AtomCache();
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
		ScopFactory.setScopDatabase(ScopFactory.VERSION_2_0_1, true);
		Sample sample = Sample.fromXML(new File(inputFilename));

		//Set up output
		PrintStream out = System.out;
		if(outputFilename != null) {
			out = new PrintStream(outputFilename);
		}
		//header
		out.print("Name\tManual\tTM\tAngle");
		for(OrderDetector detector: detectors) {
			out.print("\t");
			out.print(detector.toString());
		}
		out.println();
		
		//Recalculate orders
		for(Case c : sample.getData()) {
			CensusResult result = c.getResult();

			Atom[] ca1;
			AFPChain afpChain;
			double angle;
			try {
				ca1 = cache.getAtoms(result.getId());
				Atom[] ca2 = StructureTools.cloneAtomArray(ca1);
				afpChain = result.getAlignment().buildAfpChain(ca1, ca2);
				
				RotationAxis axis = new RotationAxis(afpChain);
				angle = axis.getAngle();
			} catch (StructureException e) {
				logger.warn("Unable to get alignment for "+result.getId(),e);
				continue;
			}


			out.print(String.format("%s\t%s\t%.3f\t%.3f",
					result.getId(),
					c.getKnownGroup(),
					result.getScoreList().getTmScore(),
					angle));
			
			for(OrderDetector detector: detectors) {
				int order = -1;
				try {
					order = detector.calculateOrder(afpChain, ca1);
				} catch (Exception e) {
					logger.warn(String.format("Skipping detector %s for %s",detector,result.getId()),e);
				}
				out.print("\t");
				out.print(order);
			}
			out.println();
		}

		// Finish up
		if(out != System.out) {
			out.close();
		}
	}

	private AtomCache cache;

	public BenchmarkOrderDetector(AtomCache cache) {
		this.cache = cache;
	}

	/**
	 * Recalculates the order for each result in sample based on the specified
	 * OrderDetector
	 * @param sample
	 * @param detector
	 * @param out An output stream to report results (may be null)
	 */
	private void calculateOrder(Sample sample, OrderDetector detector, PrintStream out) {
		OrderDetermination determination = new OrderDetectorDeterminationAdaptor(detector, cache);
		for (Case c : sample.getData()) {
			CensusResult result = c.getResult();
			int order = determination.getOrder(result);
			result.setGroup(new CensusSymmetryGroup("C" + order));
			if(out != null) {
				out.println(result.getId()+"\t"+order);
			}
		}
	}

	/**
	 * Calculates the accuracy of the order detection
	 * @param sample
	 * @return
	 */
	private double getOrderAccuracy(Sample sample ) {
		OrderAccuracy acc = new OrderAccuracy(sample, CensusSignificanceFactory.forCeSymmTm(), GroupComparisonFactory.exact());
		return acc.getAccuracy();
	}

	/**
	 * Get the accuracy of a specified OrderDetector against a sample
	 * @param sample the benchmark, annotated with true results
	 * @param detector the order detection method
	 * @return the accuracy (TP+TN)/N
	 */
	private double testOrderDetector(Sample sample, OrderDetector detector) {
		calculateOrder(sample,detector,null);
		return getOrderAccuracy(sample);
	}
}
