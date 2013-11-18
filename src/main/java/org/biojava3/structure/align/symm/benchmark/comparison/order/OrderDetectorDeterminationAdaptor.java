package org.biojava3.structure.align.symm.benchmark.comparison.order;

import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.align.symm.census2.Result;
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
	public int getOrder(Result result) {
		if (result.getAlignmentMapping() == null) {
			throw new IllegalArgumentException("Alignment mapping needed to use adaptor");
		}
		try {
			Atom[] ca = cache.getAtoms(result.getScopId());
			AFPChain afpChain = result.getAlignmentMapping().buildAfpChain(ca, StructureTools.cloneCAArray(ca));
			return detector.calculateOrder(afpChain, ca);
		} catch (StructureException e) {
			logger.error(e);
			return -1;
		} catch (IOException e) {
			logger.error(e);
			return -1;
		} catch (OrderDetectionFailedException e) {
			logger.error(e);
			return -1;
		}
	}

}
