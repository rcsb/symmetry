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
package org.biojava3.structure.align.symm.benchmark.set;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * Finds duplicates in a benchmark set (list of SCOP Ids).
 * @author dmyerstu
 */
public class DuplicateFinder {

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: DuplicateFinder input-file scop-category");
			return;
		}
		final File file = new File(args[0]);
		DuplicateFinder finder = new DuplicateFinder(ScopCategory.fromString(args[1]));
		List<String> duplicates = finder.findDuplicates(file);
		System.out.println("\n\n" + duplicates.size() + ":");
		for (String duplicate : duplicates) {
			System.out.println(duplicate);
		}
	}

	private static List<String> readNames(File file) throws IOException {
		List<String> names = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = "";
		while ((line = br.readLine()) != null) {
			if (!line.trim().equals("")) names.add(line);
		}
		br.close();
		Collections.reverse(names);
		return names;
	}

	private ScopCategory category;

	public DuplicateFinder(ScopCategory category) {
		this.category = category;
	}

	public List<String> findDuplicates(File file) throws IOException {
		return findDuplicates(readNames(file));
	}

	public List<String> findFilteredSet(File file) throws IOException {
		return findFilteredSet(readNames(file));
	}

	private List<String> findDuplicates(List<String> all) {
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75B);
		List<String> duplicates = new ArrayList<String>(all.size());
		for (int i = 0; i < all.size(); i++) {
			ScopDomain a = scop.getDomainByScopID(all.get(i));
			if (a == null) {
				System.err.println("Couldn't find " + all.get(i));
				continue;
			}
			for (int j = 0; j < i; j++) {
				ScopDomain b = scop.getDomainByScopID(all.get(j));
				if (b == null) {
					System.err.println("Couldn't find " + all.get(j));
					continue;
				}
				if (sunIdOfCategory(a) == sunIdOfCategory(b)) {
					System.out.println("(" + i + "," + j + ") " + a.getClassificationId() + "==" + b.getClassificationId());//d1l8aa3
					duplicates.add(all.get(j));
				}
			}
		}
		return duplicates;
	}

	private List<String> findFilteredSet(List<String> all) {
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75B);
		List<String> filtered = new ArrayList<String>(all.size());
		for (int i = 0; i < all.size(); i++) {
			ScopDomain a = scop.getDomainByScopID(all.get(i));
			boolean isDuplicate = false;
			for (int j = 0; j < j; j++) {
				ScopDomain b = scop.getDomainByScopID(all.get(j));
				if (sunIdOfCategory(a) == sunIdOfCategory(b)) {
					isDuplicate = true;
					break;
				}
			}
			if (!isDuplicate) {
				filtered.add(all.get(i));
			}
		}
		return filtered;
	}

	private int sunIdOfCategory(ScopDomain domain) {
		switch (category) {
		case Class:
			return domain.getClassId();
		case Fold:
			return domain.getFoldId();
		case Superfamily:
			return domain.getSuperfamilyId();
		case Family:
			return domain.getFamilyId();
		case Domain:
			return domain.getDomainId();
		default:
			throw new IllegalArgumentException("Invalid SCOP category " + category.name() + " for " + domain.getScopId() + " (" + domain.getClassificationId() + ")");
		}
	}

}
