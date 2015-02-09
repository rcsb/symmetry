package org.biojava.nbio.benchmark;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.align.client.PdbPair;

public class FisherBenchmark implements Benchmark {

	private static final String name = "Fisher, Elofsson, Rice, Eisenberg";

	List<PdbPair> pairs;
	public FisherBenchmark(){
		try {
			pairs = readFile();
		} catch (Exception e){
			e.printStackTrace();
		}
	}



	@Override
	public String getName() {
		return name;
	}

	@Override
	public List<PdbPair> getPairs() {
		return pairs;
	}


	private List<PdbPair> readFile() throws IOException {

		List<PdbPair> pairs = new ArrayList<PdbPair>();

		InputStream inStream = FisherBenchmark.class.getResourceAsStream(String.format("/fisherBenchmark.txt"));

		BufferedReader br = new BufferedReader(new InputStreamReader(inStream));
		String strLine;
		//Read File Line By Line
		while ((strLine = br.readLine()) != null)   {
			// Print the content on the console
			//System.out.println (strLine);

			String[] spl = strLine.split("\t");

			if ( spl.length != 2) {
				System.err.println("problem reading line: " + strLine);
				continue;
			}
			String name1 = spl[0];
			String name2 = spl[1];

			PdbPair pair = new PdbPair(name1, name2);
			pairs.add(pair);
		}

		return pairs;

	}


}
