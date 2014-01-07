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

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

/**
 * Determines whether a result is significant.
 * @author dmyersturnbull
 * @deprecated
 */
@Deprecated
public interface Significance {
	boolean isPossiblySignificant(AFPChain afpChain);

	boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain);
	
	boolean isSignificant(Result result);
}
