package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * Searches SCOP descriptions for particular words.
 * @author dmyerstu
 * TODO: Make work with the longer descriptions on the SCOP website
 */
public class ScopNameMiner {

	private static final Logger logger = LogManager.getLogger(ScopNameMiner.class.getName());

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length < 2) {
			System.err.println("Usage: " + ScopNameMiner.class.getSimpleName() + " orders-file-name word-1 [word-2 ...]");
			return;
		}
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75A);
		File file = new File(args[0]);
		String[] words = new String[args.length-1];
		for (int i = 0; i < args.length-1; i++) {
			words[i] = args[i+1];
		}
		ScopNameMiner miner = new ScopNameMiner(file, words);
		System.out.println(miner.printList());
	}
	
	private static String NEWLINE = "\n";
	
	private String[] words;
	private List<ScopCategory> categories = new ArrayList<ScopCategory>();
	private List<String> descriptions = new ArrayList<String>();
	private List<String> scopIds;

	public ScopNameMiner(File sampleFile, String... words) throws IOException {
		this(SampleBuilder.getNames(sampleFile), words);
	}
	
	public ScopNameMiner(List<String> scopIds, String... words) {
		this.scopIds = scopIds;
		this.words = words;
		main: for (String scopId : scopIds) {
			ScopCategory[] categories = new ScopCategory[] {ScopCategory.Fold};
			for (ScopCategory category : categories) {
				String description = getDescription(scopId, category).getDescription();
				logger.debug("Checking " + description + " of category " + category);
				for (String word : words) {
					if (description.contains(word)) {
						this.categories.add(category);
						this.descriptions.add(description);
						continue main;
					}
				}
			}
			this.categories.add(null);
			this.descriptions.add(null);
		}
	}
	
	public String printList() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < scopIds.size(); i++) {
			if (categories.get(i) == null) continue;
			sb.append(scopIds.get(i) + "\t" + categories.get(i) + "\t" + descriptions.get(i) + NEWLINE);
		}
		return sb.toString();
	}

	@SuppressWarnings("incomplete-switch")
	private static ScopDescription getDescription(String scopId, ScopCategory category) {
		ScopDatabase scop = ScopFactory.getSCOP();
		ScopDomain domain = scop.getDomainByScopID(scopId);
		if (category == ScopCategory.Domain) return scop.getScopDescriptionBySunid(domain.getSunid());
		Integer x = null;
		switch (category) {
		case Class:
			x = domain.getClassId();
		case Fold:
			x = domain.getFoldId();
		case Superfamily:
			x = domain.getSuperfamilyId();
		case Family:
			x = domain.getFamilyId();
		}
		if (x == null) return null;
		return scop.getScopDescriptionBySunid(x);
	}

}
