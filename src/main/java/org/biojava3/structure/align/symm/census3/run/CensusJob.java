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
 * Created on 2013-02-18
 */
package org.biojava3.structure.align.symm.census3.run;

import java.io.IOException;
import java.util.concurrent.Callable;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.structure.align.symm.census3.CensusAlignment;
import org.biojava3.structure.align.symm.census3.CensusAxis;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusSymmetryGroup;
import org.biojava3.structure.align.symm.census3.run.Census.AlgorithmGiver;
import org.biojava3.structure.align.symm.order.OrderDetector;
import org.biojava3.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

/**
 * One run of CE-Symm in the census.
 * 
 * @author dmyersturnbull
 */
public class CensusJob implements Callable<CensusResult> {

	private OrderDetector orderDetector = new SequenceFunctionOrderDetector();
	
	private static final Logger logger = LogManager.getLogger(CensusJob.class.getSimpleName());

	private AlgorithmGiver algorithm;
	private AfpChainCensusRestrictor significance;

	private String name;
	private Integer count;
	
	private AtomCache cache;

	private boolean storeAfpChain = false;
	private boolean recordAlignmentMapping = true; 

	private Long timeTaken;
	
	private AFPChain afpChain;

	public static CensusResult runJob(String name, int count, AlgorithmGiver algorithm, AfpChainCensusRestrictor significance, AtomCache cache) {
		return setUpJob(name, count, algorithm, significance, cache).call();
	}

	public static CensusJob setUpJob(String name, int count, AlgorithmGiver algorithm, AfpChainCensusRestrictor significance,
			AtomCache cache) {
		CensusJob job = new CensusJob(algorithm, significance);
		job.setCache(cache);
		job.setName(name);
		job.setCount(count);
		return job;
	}

	private static boolean sanityCheckPreAlign(Atom[] ca1, Atom[] ca2) {
		if (ca1 == ca2) return false;
		if (ca1[0].getGroup().getChain().getParent() == ca2[0].getGroup().getChain().getParent()) return false;
		return true;
	}

	public CensusJob(AlgorithmGiver algorithm, AfpChainCensusRestrictor significance) {
		this.algorithm = algorithm;
		this.significance = significance;
	}

	@Override
	public CensusResult call() {

		if (name == null || count == null) throw new IllegalStateException("Must set name and count first.");

		// first, get the atoms
		Atom[] ca1, ca2;
		Structure structure;
		logger.debug("Getting atoms for " + name + " (job #" + count + ")");
		try {
			if (cache == null) cache = new AtomCache();
			structure = cache.getStructure(name);
			ca1 = StructureTools.getAtomCAArray(structure);
			ca2 = StructureTools.cloneCAArray(ca1);
		} catch (Exception e) {
			logger.error("Could not create the atom arrays for " + name + ": " + e.getMessage(), e);
			return null;
		}
		logger.debug("Got " + ca1.length + " atoms (job #" + count + ")");

		// run the alignment
		StructureAlignment algorithm = this.algorithm.getAlgorithm();
		AFPChain afpChain = null;
		logger.debug("Running CE-Symm (job #" + count + ")");
		try {
			afpChain = findSymmetry(name, algorithm, ca1, ca2);
		} catch (Exception e) {
			logger.error("Failed running CE-Symm on " + name + ": " + e.getMessage(), e);
			return convertResult(null, null, null, null, name, null);
		}
		if (afpChain == null || afpChain.getOptAln() == null) {
			logger.debug("CE-Symm returned null (job #" + count + ")");
			return convertResult(null, null, null, null, name, null);
		}

		// there are two cases in which we know there is no symmetry
		if (afpChain.getBlockNum() != 2) {
			logger.debug("CE-Symm returned a result with " + afpChain.getBlockNum() + " block(s) (job #" + count + ")");
			return convertResult(afpChain, false, null, null, name, null);
		}
		if (afpChain.getAlnLength() < 1) {
			logger.debug("CE-Symm returned an empty alignment (job #" + count + ")");
			return convertResult(afpChain, false, null, null, name, null);
		}

		if (significance.isPossiblySignificant(afpChain)) {

			logger.debug("Result is significant (job #" + count + ")");

			// initializes these to null, then try setting them separately
			Boolean isSignificant = null;
			Integer order = null;
			Float angle = null;
			Protodomain protodomain = null;

			// first try to find the protodomain
			logger.debug("Finding protodomain (job #" + count + ")");
			try {
				protodomain = Protodomain.fromSymmetryAlignment(afpChain, ca2, 1, cache);
				logger.debug("Protodomain is " + protodomain + " (job #" + count + ")");
			} catch (Exception e) {
				logger.warn("Could not create protodomain because " + e.getMessage(), e);
			}

			// now try to find the order
			logger.debug("Finding order (job #" + count + ")");
			try {
				order = orderDetector.calculateOrder(afpChain, ca1);
				logger.debug("Order is " + order + " (job #" + count + ")");
			} catch (Exception e) {
				logger.error("Failed to determine the order of symmetry on " + name + ": " + e.getMessage(), e);
			}

			// now try to find the angle
			logger.debug("Finding angle (job #" + count + ")");
			try {
				angle = (float) getAngle(afpChain, ca1, ca2);
				logger.debug("Angle is " + angle + " (job #" + count + ")");
			} catch (Exception e) {
				logger.error("Failed to determine the angle on " + name + ": " + e.getMessage(), e);
			}

			return convertResult(afpChain, isSignificant, order, protodomain, name, angle);

		} else { // trying this would take too long
			logger.debug("Result is not significant (job #" + count + ")");
			return convertResult(afpChain, false, null, null, name, null);
		}
	}

	public AFPChain getAfpChain() {
		return afpChain;
	}

	public Long getTimeTaken() {
		return timeTaken;
	}

	public void setRecordAlignmentMapping(boolean recordAlignmentMapping) {
		this.recordAlignmentMapping = recordAlignmentMapping;
	}

	public void setCache(AtomCache cache) {
		this.cache = cache;
	}

	public void setCount(int count) {
		this.count = count;
	}

	public void setName(String name) {
		this.name = name;
	}

	public void setOrderDetector(OrderDetector orderDetector) {
		this.orderDetector = orderDetector;
	}

	public void nullifyAfpChain() {
		this.afpChain = null;
	}
	
	/**
	 * Store the raw AFPChain in {@link #getAfpChain()}. Requires a lot of memory. Not recommended.
	 * 
	 * @param storeAfpChain
	 */
	public void setStoreAfpChain(boolean storeAfpChain) {
		this.storeAfpChain = storeAfpChain;
	}

	private CensusResult convertResult(AFPChain afpChain, Boolean isSymmetric, Integer order, Protodomain protodomain,
			String name, Float angle) {

		CensusResult r = new CensusResult();

		if (storeAfpChain) this.afpChain = afpChain;
		if (recordAlignmentMapping) r.setAlignment(new CensusAlignment(afpChain));

		if (protodomain != null) r.setAlignedUnit(protodomain.toString());

		r.setAlignment(new CensusAlignment(afpChain));
		try {
			if (afpChain.getAlnLength() > 0) r.setAxis(new CensusAxis(new RotationAxis(afpChain)));

		} catch (RuntimeException e) {

			logger.error("Could not get rotation axis for " + name + "(job #" + count + ")", e);

		} catch (Exception e) {

			e.printStackTrace();

			logger.error("Alignment for " + name + " is empty (job #" + count + ")", e);

			if (angle != null) { // if the axis can't be found, at least we do have the angle
				CensusAxis axis = new CensusAxis();
				axis.setAngle((float) angle);
				r.setAxis(axis);
			}
		}
		r.setGroup(new CensusSymmetryGroup("C" + order));

		return r;
	}

	private AFPChain findSymmetry(String name, StructureAlignment algorithm, Atom[] ca1, Atom[] ca2) throws StructureException, IOException {
		if (!sanityCheckPreAlign(ca1, ca2)) throw new RuntimeException("Can't align using same structure.");
		long startTime = System.currentTimeMillis();
		AFPChain afpChain = algorithm.align(ca1, ca2);
		long endTime = System.currentTimeMillis();
		timeTaken = endTime - startTime;
		if (afpChain == null) return null;
		afpChain.setName1(name);
		afpChain.setName2(name);
		double realTmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
		afpChain.setTMScore(realTmScore);
		return afpChain;
	}

	/**
	 * Returns the <em>magnitude</em> of the angle between the first and second blocks of {@code afpChain}, measured in
	 * degrees. This is always a positive value (unsigned).
	 * 
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	private double getAngle(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {
		Matrix rotation = afpChain.getBlockRotationMatrix()[0];
		return Math.acos(rotation.trace() - 1) * 180 / Math.PI;
	}

}
