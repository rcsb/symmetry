package org.biojava3.structure.quaternary.misc;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SimpleCsvReader {
	private List<String> header = Collections.emptyList();
	private List<List<String>> csv = new ArrayList<List<String>>();
	
	public void readFile(String fileName) throws IOException {
		header = Collections.emptyList();
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		String line = null;
		line = reader.readLine();
		if (line != null) {
			header = Arrays.asList(line.split(","));
		}
		while ((line = reader.readLine()) != null) {
			String[] tokens = line.split(",");
			csv.add(Arrays.asList(tokens));
		}
		reader.close();
	}
	
	public List<String> getRow(int rowIndex) {
		return csv.get(rowIndex);
	}
	
	public List<String> getColumn(String columName) {
		List<String> column = new ArrayList<String>(csv.size());
		int index = header.indexOf(columName);
		if (index >= 0) {
			for (List<String> row: csv) {
				column.add(row.get(index));
			}
		}
		return column;
	}
}
