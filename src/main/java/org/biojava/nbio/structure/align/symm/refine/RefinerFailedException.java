package org.biojava.nbio.structure.align.symm.refine;

/**
 * Refine alignment failed.
 * @author lafita
 */
public class RefinerFailedException extends Exception {
	private static final long serialVersionUID = -7040421412578699838L;

	public RefinerFailedException() {
		super();
	}

	public RefinerFailedException(String message, Throwable cause) {
		super(message, cause);
	}

	public RefinerFailedException(String message) {
		super(message);
	}

	public RefinerFailedException(Throwable cause) {
		super(cause);
	}

}
