package org.biojava.nbio.structure.align.symm.census3.analysis.mespeus;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.biojava.nbio.structure.align.symm.census3.CensusResult;
import org.biojava.nbio.structure.align.symm.census3.CensusResultList;
import org.biojava.nbio.structure.align.symm.census3.CensusSignificance;
import org.biojava.nbio.structure.align.symm.census3.CensusSignificanceFactory;

/**
 * Quantifies symmetry among metalloproteins in <a href="http://mespeus.bch.ed.ac.uk/MESPEUS">MESPEUS</a>.
 * @author dmyersturnbull
 */
public class SymmetryInMespeus {


	private CensusSignificance significance;
	private List<MespeusEntry> entries;

	public List<MespeusEntry> getEntries() {
		return entries;
	}

	public SymmetryInMespeus(File mespeus, CensusSignificance significance) throws IOException {
		if (significance == null) throw new IllegalArgumentException("Significance can't be null");
		this.significance = significance;
		entries = new ArrayList<MespeusEntry>();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(mespeus));
			br.readLine(); // skip title row
			String line = "";
			while ((line = br.readLine()) != null) {
				entries.add(MespeusEntry.parse(line));
			}
		} finally {
			if (br != null) br.close();
		}
	}

	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: " + SymmetryInMespeus.class.getSimpleName() + " census-file.xml mespeus-database-file.tsv");
			return;
		}
		CensusSignificance significance = CensusSignificanceFactory.forCeSymmOrd();
		SymmetryInMespeus mespeus = new SymmetryInMespeus(new File(args[1]), significance);
		CensusResultList census = CensusResultList.fromXML(new File(args[0]));
		DescriptiveStatistics stats = mespeus.correlate(census, MespeusEntryMatcherFactory.everything());
//		DescriptiveStatistics stats = mespeus.correlate(census, MespeusEntryMatcher.metalName(new String[] {"FE", "FE1", "FE2", "FE3", "FE4"}));
		System.out.println(stats);
		DescriptiveStatistics symmStats = mespeus.correlateWithSymmetryOrder(census, MespeusEntryMatcherFactory.everything());
		System.out.println(symmStats);
	}

	private DescriptiveStatistics correlateWithSymmetryOrder(CensusResultList census, MespeusEntryMatcher everything) {
		Map<String,Integer> orders = new HashMap<String,Integer>();
		for (CensusResult result : census.getEntries()) {
			if (result.getAlignedUnit() == null) continue;
			int x = result.getAlignedUnit().indexOf('.');
			if (x == -1) continue;
			String pdbId = result.getAlignedUnit().substring(0, x);
			if (!orders.containsKey(pdbId)) {
				if (significance.isSignificant(result)) {
					orders.put(pdbId, result.getOrder());
				} else {
					orders.put(pdbId, 1);
				}
			}
		}

		DescriptiveStatistics stats = new DescriptiveStatistics();
		for (MespeusEntry entry : entries) {
			if (everything.matches(entry)) {
				Integer order = orders.get(entry.getPdbId());
				if (order != null) {
					MespeusEntryMatcher matcher = MespeusEntryMatcherFactory.coordinationEquals(order);
					stats.addValue(matcher.matches(entry)? 1 : 0);
				}
			}
		}
		return stats;
	}

	public DescriptiveStatistics correlate(CensusResultList census, MespeusEntryMatcher matcher) {
		Map<String,Integer> orders = new HashMap<String,Integer>();
		for (CensusResult result : census.getEntries()) {
			if (result.getAlignedUnit() == null) continue;
			int x = result.getAlignedUnit().indexOf('.');
			if (x == -1) continue;
			String pdbId = result.getAlignedUnit().substring(0, x);
			if (!orders.containsKey(pdbId)) {
				if (significance.isSignificant(result)) {
					orders.put(pdbId, result.getOrder());
				} else {
					orders.put(pdbId, 1);
				}
			}
		}

		DescriptiveStatistics stats = new DescriptiveStatistics();
		for (MespeusEntry entry : entries) {
			if (matcher.matches(entry)) {
				Integer order = orders.get(entry.getPdbId());
				if (order != null) {
					stats.addValue(order>1? 1 : 0);
				}
			}
		}
		return stats;
	}

}