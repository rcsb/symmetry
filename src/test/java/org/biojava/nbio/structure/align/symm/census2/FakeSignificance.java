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
 * Created on 2013-03-22
 *
 */
package org.biojava.nbio.structure.align.symm.census2;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.census2.Result;
import org.biojava.nbio.structure.align.symm.census2.Significance;
import org.biojava.nbio.structure.align.symm.protodomain.Protodomain;

/**
 * A mock {@link Significance} object for use by {@link SignificanceFactoryTest}. This should be a seperate file to make calling with reflection easier.
 * @author dmyerstu
 */
public class FakeSignificance implements Significance {

	public FakeSignificance() {
		
	}

	public FakeSignificance(String string) {
		
	}
	
	@Override
	public boolean isPossiblySignificant(AFPChain afpChain) {
		return false;
	}

	@Override
	public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
		return false;
	}

	@Override
	public boolean isSignificant(Result result) {
		return false;
	}

}
