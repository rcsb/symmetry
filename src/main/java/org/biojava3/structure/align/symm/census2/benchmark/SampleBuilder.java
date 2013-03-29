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
 * Created on 2013-03-03
 *
 */
package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;

/**
 * Converts a {@link Results} to a {@link Sample}.
 * @author dmyerstu
 */
public class SampleBuilder {

	public static void main(String[] args) throws IOException {
		File input = new File(args[0]);
		File output = new File(args[1]);
		File ordersFile = new File(args[2]);
		buildSample(input, output, getOrders(ordersFile));
	}
	
	
	public static Map<String,KnownInfo> getOrders(File knownInfoFile) throws IOException {
		Map<String,KnownInfo> map = new HashMap<String,KnownInfo>();
		BufferedReader br = new BufferedReader(new FileReader(knownInfoFile));
		String line = "";
		while ((line = br.readLine()) != null) {
			String[] parts = line.split("\t");
			KnownInfo info = new KnownInfo(parts[2], Integer.parseInt(parts[1]));
			map.put(parts[0], info);
		}
		br.close();
		return map;
	}
	
	public static void buildSample(File input, File output, Map<String,KnownInfo> knownInfos) throws IOException {
		Results results = Results.fromXML(input);
		Sample sample = new Sample();
		for (Result result : results.getData()) {
			if (result == null) continue;
			Case c = new Case();
			c.setResult(result);
			c.setKnownInfo(knownInfos.get(result.getScopId()));
			sample.add(c);
		}
		BufferedWriter br = new BufferedWriter(new FileWriter(output));
		br.write(sample.toXML());
		br.close();
	}


	public static Map<String, KnownInfo> getOrders(String file) throws IOException {
		return getOrders(new File(file));
	}


	public static void buildSample(File input, File output, File ordersFile) throws IOException {
		buildSample(input, output, getOrders(ordersFile.getPath()));
	}

}
