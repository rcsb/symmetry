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
 * Created on 2013-05-23
 *
 */
package org.biojava3.structure.quaternary.core;

/**
 * Holds the results of quaternary symmetry perception.
 * 
 * @author Peter Rose
 *
 */
public class QuatSymmetryResults {
	private Subunits subunits = null;
	private RotationGroup rotationGroup = null;
	private String method = null;
	private double sequenceIdentityThreshold = 0;
	private boolean local = false;
	private boolean preferredResult = false;
	
	public QuatSymmetryResults(Subunits subunits, RotationGroup rotationGroup, String method) {
		this.subunits = subunits;
		this.rotationGroup = rotationGroup;
		this.method = method;
	}
	
	/**
	 * Returns protein subunit information that was used to determine symmetry information
	 * 
	 * @return
	 */
	public Subunits getSubunits() {
		return subunits;
	}
	
	/**
	 * Returns rotation group (point group) information representing rotational quaternary symmetry,
	 * see http://en.wikipedia.org/wiki/Rotation_group_SO(3)
	 * 
	 * @return rotation group
	 */
	public RotationGroup getRotationGroup() {
		return rotationGroup;
	}

	/**
	 * Returns name of method used for symmetry perception.
	 * 
	 * @return method
	 */
	public String getMethod() {
		return method;
	}

	public double getSequenceIdentityThreshold() {
		return sequenceIdentityThreshold;
	}

	public void setSequenceIdentityThreshold(double sequenceIdentityThreshold) {
		this.sequenceIdentityThreshold = sequenceIdentityThreshold;
	}

	/**
	 * Return true 
	 * @return
	 */
	public boolean isLocal() {
		return local;
	}

	public void setLocal(boolean local) {
		this.local = local;
	}

	public boolean isPreferredResult() {
		return preferredResult;
	}

	public void setPreferredResult(boolean preferredResult) {
		this.preferredResult = preferredResult;
	}
}
