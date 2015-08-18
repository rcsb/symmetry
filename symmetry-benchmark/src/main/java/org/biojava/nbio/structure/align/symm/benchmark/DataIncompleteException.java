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
 * Created on 2013-02-22
 *
 */
package org.biojava.nbio.structure.align.symm.benchmark;

/**
 * Data required for the benchmark was not found.
 * @author dmyerstu
 */
public class DataIncompleteException extends RuntimeException {

	private static final long serialVersionUID = 7865520194833163469L;
	private String scopId;
	
	public String getScopId() {
		return scopId;
	}

	public DataIncompleteException(String scopId) {
		super();
		this.scopId = scopId;
	}

	public DataIncompleteException(String scopId, String message) {
		super(message);
		this.scopId = scopId;
	}

	public DataIncompleteException(String scopId, Throwable cause) {
		super(cause);
		this.scopId = scopId;
	}

	public DataIncompleteException(String scopId, String message, Throwable cause) {
		super(message, cause);
		this.scopId = scopId;
	}

	
}
