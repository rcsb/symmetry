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
package org.biojava3.structure.align.symm.census2.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava3.structure.align.symm.census2.Results;

/**
 * A utility that combines multiple census files into one census file.
 * @author dmyerstu
 */
public class CensusCombiner {

	private Results results;

	public CensusCombiner(File[] files) {
		try {
			results = Results.fromXML(files);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public static void main(String[] args) {
		File combined = new File(args[0]);
		File[] files = new File[args.length-1];
		for (int i = 1; i < args.length; i++) {
			files[i-1] = new File(args[i]);
		}
		CensusCombiner combiner = new CensusCombiner(files);
		combiner.print(combined);
	}

	public void print(File combined) {
		PrintWriter out = null;
		try {
			out = new PrintWriter(new BufferedWriter(new FileWriter(combined)));
			String xml;
			xml = results.toXML();
			out.print(xml);
			out.flush();
			out.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			out.close();
		}
	}

}
