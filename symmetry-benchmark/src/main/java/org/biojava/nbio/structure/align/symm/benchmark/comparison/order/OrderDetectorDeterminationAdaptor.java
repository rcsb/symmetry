package org.biojava.nbio.structure.align.symm.benchmark.comparison.order;



import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.census3.CensusResult;
import org.biojava.nbio.structure.align.symm.order.OrderDetector;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * An adaptor for symmetry-benchmark order analysis code ({@link OrderDetermination}) that uses order-detection code in the symmetry project ({@link OrderDetector}).
 * @author dmyersturnbull
 */
public class OrderDetectorDeterminationAdaptor implements OrderDetermination {

	private static final Logger logger = LoggerFactory.getLogger(OrderDetectorDeterminationAdaptor.class);

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
			if(result.getId() == null) {
				throw new NullPointerException("Null ID");
			}
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
