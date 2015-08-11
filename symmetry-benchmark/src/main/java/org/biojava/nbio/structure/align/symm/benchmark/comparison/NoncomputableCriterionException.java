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
package org.biojava.nbio.structure.align.symm.benchmark.comparison;

/**
 * The {@link Criterion} could not determine the quality of the {@link Result} because the result was missing critical data.
 * @author dmyerstu
 */
public class NoncomputableCriterionException extends Exception {

	private static final long serialVersionUID = 5469531500394989349L;

	public NoncomputableCriterionException() {
		super();
	}

	public NoncomputableCriterionException(String message, Throwable cause) {
		super(message, cause);
	}

	public NoncomputableCriterionException(String message) {
		super(message);
	}

	public NoncomputableCriterionException(Throwable cause) {
		super(cause);
	}

}
