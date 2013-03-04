/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 2013-02-18
 *
 */
package org.biojava3.structure.align.symm.census2;

import java.io.IOException;
import java.util.concurrent.Callable;

import org.apache.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.census2.Census.AlgorithmGiver;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

/**
 * One run of CE-Symm in the census.
 * @author dmyerstu
 */
public class CensusJob implements Callable<Result> {

	static final Logger logger = Census.logger;

	private AlgorithmGiver algorithm;

	private AtomCache cache;
	private Integer count;

	private ScopDomain domain;
	private Significance significance;
	private ScopDescription superfamily;

	private static boolean sanityCheckPreAlign(Atom[] ca1, Atom[] ca2) {
		if (ca1 == ca2) return false;
		if (ca1[0].getGroup().getChain().getParent() == ca2[0].getGroup().getChain().getParent()) return false;
		return true;
	}

	public CensusJob(AtomCache cache, AlgorithmGiver algorithm, Significance significance) {
		this.algorithm = algorithm;
		this.cache = cache;
		this.significance = significance;
	}

	@Override
	public Result call() {
		if (domain == null || superfamily == null || count == null) throw new IllegalStateException (
				"Must set domain, superfamily, and count first.");
		String name = domain.getScopId();
		Atom[] ca1, ca2;
		logger.debug("Getting atoms for " + name + " (job #" + count + ")");
		try {
			ca1 = cache.getAtoms(name);
			ca2 = StructureTools.cloneCAArray(ca1);
		} catch (Exception e) {
			logger.error("Could not create the atom arrays for " + name + ": " + e.getMessage(), e);
			return null;
		}
		logger.debug("Got " + ca1.length + " atoms (job #" + count + ")");
		AFPChain afpChain;
		logger.debug("Running CE-Symm (job #" + count + ")");
		try {
			afpChain = findSymmetry(name, ca1, ca2);
		} catch (Exception e) {
			logger.error("Failed running CE-Symm on " + name + ": " + e.getMessage(), e);
			return convertResult(null, null, superfamily, name, null, null, domain);
		}
		if (afpChain == null || afpChain.getOptAln() == null) {
			logger.debug("CE-Symm returned null (job #" + count + ")");
			return convertResult(null, null, superfamily, name, null, null, domain);
		}
		if (afpChain.getBlockNum() != 2) {
			logger.debug("CE-Symm returned a result with " + afpChain.getBlockNum() + " block(s) (job #" + count + ")");
			return convertResult(afpChain, false, superfamily, name, null, null, domain);
		}

		if (significance.isPossiblySignificant(afpChain)) {
			logger.debug("Result is significant (job #" + count + ")");
			logger.debug("Finding protodomain (job #" + count + ")");
			Protodomain protodomain;
			try {
				protodomain = Protodomain.fromSymmetryAlignment(afpChain, ca2, 1, cache);
				logger.debug("Protodomain is " + protodomain + " (job #" + count + ")");
			} catch (Exception e) {
				logger.warn("Could not create protodomain because " + e.getMessage(), e);
				return convertResult(afpChain, false, superfamily, name, null, null, domain);
			}
			logger.debug("Finding order (job #" + count + ")");
			int order;
			try {
				order = CeSymm.getSymmetryOrder(afpChain);
				logger.debug("Order is " + order + " (job #" + count + ")");
			} catch (Exception e) {
				logger.error("Failed to determine the order of symmetry on " + name + ": " + e.getMessage(), e);
				return convertResult(afpChain, false, superfamily, name, null, null, domain);
			}
			logger.debug("Finding angle (job #" + count + ")");
			double angle;
			try {
				angle = getAngle(afpChain, ca1, ca2);
				logger.debug("Angle is " + angle + " (job #" + count + ")");
			} catch (Exception e) {
				logger.error("Failed to determine the angle on " + name + ": " + e.getMessage(), e);
				return convertResult(null, false, superfamily, name, order, null, domain);
			}
			boolean isSymmetric = significance.isSignificant(protodomain, order, angle, afpChain);
			return convertResult(afpChain, isSymmetric, superfamily, name, order, protodomain.toString(),
					domain);
		} else {
			logger.debug("Result is not significant (job #" + count + ")");
			return convertResult(afpChain, false, superfamily, name, null, null, domain);
		}
	}

	/**
	 * Returns the <em>magnitude</em> of the angle between the first and second blocks of {@code afpChain}, measured in degrees. This is always a positive value (unsigned).
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	private double getAngle(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {
		Matrix rotation = afpChain.getBlockRotationMatrix()[0];
		return Math.acos(rotation.trace() - 1) * 180/Math.PI;
	}

	public void setCount(int count) {
		this.count = count;
	}

	public void setDomain(ScopDomain domain) {
		this.domain = domain;
	}

	public void setSuperfamily(ScopDescription superfamily) {
		this.superfamily = superfamily;
	}

	private Result convertResult(AFPChain afpChain, Boolean isSymmetric, ScopDescription superfamily, String scopId, Integer order, String protodomain, ScopDomain domain) {

		final String description = superfamily.getDescription();

		Result r = new Result();

		r.setRank(count);
		r.setScopId(scopId);
		r.setClassification(superfamily.getClassificationId());
		r.setDescription(description);
		r.setProtodomain(protodomain);
		r.setSunId(domain.getSunid());

		r.setAlignment(new Alignment(afpChain));
		r.setIsSignificant(isSymmetric);
		try {
			if (afpChain != null) r.setAxis(new Axis(new RotationAxis(afpChain)));
		} catch (RuntimeException e) {
			logger.error("Could not get rotation axis for " + scopId, e);
		}
		r.setOrder(order);

		return r;
	}

	private AFPChain findSymmetry(String name, Atom[] ca1, Atom[] ca2) throws StructureException, IOException {
		if (!sanityCheckPreAlign(ca1, ca2)) throw new RuntimeException("Can't align using same structure.");
		AFPChain afpChain = algorithm.getAlgorithm().align(ca1, ca2);
		if (afpChain == null) return null;
		afpChain.setName1(name);
		afpChain.setName2(name);
		double realTmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
		afpChain.setTMScore(realTmScore);
		return afpChain;
	}

}
