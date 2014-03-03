package org.biojava3.structure.align.symm.benchmark.comparison.order;

import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava3.structure.align.symm.order.OrderDetector;

/**
 * An adaptor for symmetry-benchmark order analysis code ({@link OrderDetermination}) that uses order-detection code in the symmetry project ({@link OrderDetector}).
 * @author dmyersturnbull
 */
public class OrderDetectorDeterminationAdaptor implements OrderDetermination {

	private static final Logger logger = LogManager.getLogger(OrderDetectorDeterminationAdaptor.class.getName());

	private OrderDetector detector;
	private AtomCache cache;

	public OrderDetectorDeterminationAdaptor(OrderDetector detector, AtomCache cache) {
		this.detector = detector;
		this.cache = cache;
	}

	@Override
	public int getOrder(CensusResult result) {
		if (result.getAlignment() == null) {
			throw new IllegalArgumentException("Alignment mapping needed to use adaptor");
		}
		try {
			Atom[] ca1 = cache.getAtoms(result.getId());
			Atom[] ca2 = cache.getAtoms(result.getId());
			AFPChain afpChain = result.getAlignment().buildAfpChain(ca1, ca2);
			return detector.calculateOrder(afpChain, ca1);
		} catch (Exception e) {
			logger.error("Failed to get AFPChain from " + result.getId(), e);
			return -1;
		}
	}

}
