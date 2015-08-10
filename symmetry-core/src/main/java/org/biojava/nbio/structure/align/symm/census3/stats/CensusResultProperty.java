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
 * Created on 2013-03-10
 *
 */
package org.biojava.nbio.structure.align.symm.census3.stats;

import org.biojava.nbio.structure.align.symm.census3.CensusResult;


/**
 * A property about a {@link CensusResult}. Usually from {@link CensusResult#getScoreList()}.
 * @author dmyersturnbull
 */
public interface CensusResultProperty<T extends Number> {

	T getProperty(CensusResult result) throws PropertyUndefinedException;
	
	String getName();
	
}
