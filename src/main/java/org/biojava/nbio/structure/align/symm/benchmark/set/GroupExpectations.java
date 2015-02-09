package org.biojava.nbio.structure.align.symm.benchmark.set;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.align.symm.benchmark.KnownInfo;
import org.biojava.nbio.structure.align.symm.benchmark.SampleBuilder;

/**
 * Checks the benchmark set against per-fold or per-superfamily expectations for symmetry.
 * 
 * @author dmyerstu
 */
public class GroupExpectations {

	private static final String NEWLINE = "\n";

	private Map<String, KnownInfo> actual;
	private LinkedHashMap<String, String> expected;
	private LinkedHashMap<String, String> wrong = new LinkedHashMap<String, String>();

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length > 2) {
			System.err.println("Usage: " + GroupExpectations.class.getSimpleName() + " [actual-file] [expected-file]");
		}
		File actual = new File("src/main/resources/domain_symm_benchmark.tsv");
		File expected = new File("src/main/resources/expected_groups.tsv");
		if (args.length > 0) actual = new File(args[0]);
		if (args.length > 1) expected = new File(args[1]);
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75A);
		GroupExpectations expectations = new GroupExpectations(actual, expected);
		System.out.println(expectations);
		System.out.println("=======");
		System.out.println(expectations.toLongList());
	}

	public GroupExpectations(File actual, File expected) throws IOException {
		this.actual = SampleBuilder.getOrders(actual);
		this.expected = new LinkedHashMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(expected));
		String line = "";
		while ((line = br.readLine()) != null) {
			String[] parts = line.split("\t");
			this.expected.put(parts[0], parts[1]);
		}
		br.close();
		init();
	}

	public GroupExpectations(LinkedHashMap<String, KnownInfo> actual, LinkedHashMap<String, String> expected) {
		this.actual = actual;
		this.expected = expected;
		init();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String, String> entry : wrong.entrySet()) {
			sb.append(entry.getKey() + "\t"
					+ ScopFactory.getSCOP().getDomainByScopID(entry.getKey()).getClassificationId() + "\t"
					+ entry.getValue() + NEWLINE);
		}
		return sb.toString();
	}

	public String toLongList() {
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String, KnownInfo> entry : actual.entrySet()) {
			if (wrong.containsKey(entry.getKey())) {
				sb.append("\tsuspect");
			}
			sb.append(NEWLINE);
		}
		return sb.toString();
	}
	
	private void init() {
		ScopDatabase scop = ScopFactory.getSCOP();
		for (Map.Entry<String, KnownInfo> entry : actual.entrySet()) {
			String scopId = entry.getKey();
			String group = entry.getValue().getGroup();
			ScopDomain domain = scop.getDomainByScopID(scopId);
			String cf = expected.get(scop.getScopDescriptionBySunid(domain.getFoldId()).getClassificationId());
			if (cf != null && !group.equals(cf)) wrong.put(scopId, group);
			String sf = expected.get(scop.getScopDescriptionBySunid(domain.getFoldId()).getClassificationId());
			if (sf != null && !group.equals(cf)) wrong.put(scopId, group);
			String fa = expected.get(scop.getScopDescriptionBySunid(domain.getFoldId()).getClassificationId());
			if (fa != null && !group.equals(cf)) wrong.put(scopId, group);
		}
	}

}
