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
 * Created on 2013-02-20
 *
 */
package org.biojava3.structure.align.symm.census3.utils;

import java.io.File;
import java.io.IOException;

import org.biojava3.core.sequence.io.util.IOUtils;
import org.biojava3.structure.align.symm.census3.CensusResultList;

/**
 * A utility that combines multiple census files into one census file.
 * @author dmyersturnbull
 */
public class CensusCombiner {

	private CensusResultList results;

	public CensusCombiner(File[] files) {
		try {
			results = CensusResultList.fromXML(files);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public static void main(String[] args) {
		if (args.length < 2) {
			System.err.println("Usage: " + CensusCombiner.class.getSimpleName() + " combined-output-file.xml input-file-1.xml [input-file-2.xml ...]");
			return;
		}
		File combined = new File(args[0]);
		File[] files = new File[args.length-1];
		for (int i = 1; i < args.length; i++) {
			files[i-1] = new File(args[i]);
		}
		CensusCombiner combiner = new CensusCombiner(files);
		combiner.print(combined);
	}

	public void print(File combined) {
		try {
			Census2Adaptor.print(results.toXML(), combined);
		} catch (IOException e) {
			throw new RuntimeException("Couldn't print to file " + combined.getPath(), e);
		}
	}

}
