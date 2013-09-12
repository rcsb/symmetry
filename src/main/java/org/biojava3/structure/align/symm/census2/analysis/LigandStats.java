package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;

/**
 * Prints statistics about {@link LigandList LigandLists}.
 * @author dmyerstu
 */
public class LigandStats {

	public static void printStats(LigandList list) {
		int nMetallic = 0;
		for (StructureLigands inStruct : list.getData().values()) {
			for (Ligand ligand : inStruct.getLigands()) {
				if (ligand.isMetallic()) {
					nMetallic++;
					break;
				}
			}
		}
		System.out.println(nMetallic + " / " + list.size());
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("Usage: " + LigandFinder.class.getSimpleName() + " ligand-file.xml");
			return;
		}
		LigandList list = LigandList.fromXml(new File(args[0]));
		printStats(list);
	}

}
