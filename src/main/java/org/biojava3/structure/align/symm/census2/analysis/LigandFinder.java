package org.biojava3.structure.align.symm.census2.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomPositionMap;
import org.biojava.bio.structure.AtomPositionMap.GroupMatcher;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.census2.stats.StatUtils;

/**
 * Find ligands near the centroids of symmetric domains.
 * @author dmyerstu
 */
public class LigandFinder {

	public static int DEFAULT_RADIUS = 5;
	
	private static final Logger logger = LogManager.getLogger(LigandFinder.class.getName());

	GroupMatcher exclusionMatcher = new GroupMatcher() {
		@Override
		public boolean matches(Group group) {
			ResidueType type = group.getChemComp().getResidueType();
			return group.hasAtom(StructureTools.caAtomName) || AtomPositionMap.AMINO_ACID_NAMES.contains(group.getPDBName()) || type == ResidueType.lPeptideLinking || type == ResidueType.glycine || type == ResidueType.lPeptideAminoTerminus || type == ResidueType.lPeptideCarboxyTerminus || type == ResidueType.dPeptideLinking || type == ResidueType.dPeptideAminoTerminus || type == ResidueType.dPeptideCarboxyTerminus;
		}
	};
	private int radius = DEFAULT_RADIUS;
	private File output;
	private int printFrequency = 20;

	public void setPrintFrequency(int printFrequency) {
		this.printFrequency = printFrequency;
	}

	private Map<String,String> formulas = new LinkedHashMap<String,String>();

	public LigandFinder(int radius) {
		this.radius = radius;
	}

	public LigandFinder(int radius, GroupMatcher exclusionMatcher) {
		this.radius = radius;
		this.exclusionMatcher = exclusionMatcher;
	}

	public void setExclusionMatcher(GroupMatcher exclusionMatcher) {
		this.exclusionMatcher = exclusionMatcher;
	}

	public void setRadius(int radius) {
		this.radius = radius;
	}

	public void setOutput(File output) {
		this.output = output;
	}

	private String getEmpericalFormula(List<Atom> atoms) {
		SortedMap<String,Integer> letters = new TreeMap<String,Integer>();
		boolean hasMetal = false;
		for (Atom atom : atoms) {
			String name = atom.getElement().name();
			if (atom.getElement().getElementType().isMetal()) hasMetal = true;
			StatUtils.plus(letters, name);
		}
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String,Integer> entry : letters.entrySet()) {
			sb.append(entry.getKey());
			if (entry.getValue() > 1) sb.append(entry.getValue());
		}
		return hasMetal? "*" + sb.toString() : sb.toString();
	}

	public void find(Results census) {

		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
		Significance sig = SignificanceFactory.rotationallySymmetricSmart();
		AtomCache cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);

		int i = 0;
		for (Result result : census.getData()) {

			String scopId = result.getScopId();
			try {

				ScopDomain domain = scop.getDomainByScopID(scopId);
				if (domain == null) {
					logger.error(result.getScopId() + " is null");
					continue;
				}
				if (!sig.isSignificant(result)) {
					continue;
				}

				// we want to get all groups in the center that are NOT amino acids or water atoms
				Structure structure = cache.getStructureForDomain(scopId, scop);
				Atom[] ca = StructureTools.getAtomCAArray(structure);
				Atom centroid = Calc.getCentroid(ca);
				AtomPositionMap atomPositions = new AtomPositionMap(ca, exclusionMatcher);
				Set<ResidueNumber> residues = atomPositions.getNavMap().keySet();
				Set<Group> ligands = StructureTools.getGroupsWithinShell(structure, centroid, residues, radius, false);

				if (!ligands.isEmpty()) {
					//					System.out.println(scopId);
					for (Group group : ligands) {
						String formula = getEmpericalFormula(group.getAtoms());
						if (formulas.containsKey(scopId)) {
							formulas.put(scopId, formulas.get(scopId) + ", " + formula);
						} else {
							formulas.put(scopId, formula);
						}
						//						System.out.println("\t" + group.getType());
						//						System.out.println("\t" + getEmpericalFormula(group.getAtoms()));
						//						System.out.print("\t");
						//						for (Atom atom : group.getAtoms()) {
						//							System.out.print(" " + atom.getElement());
						//							boolean isMetal = atom.getElement().getElementType().isMetal();
						//							if (isMetal) System.out.print("*");
						//						}
						//						System.out.println();
						if (i % printFrequency == 0 && i > 0 && output != null) {
							logger.info("Printing to " + output.getPath());
							writeToFile();
						}
					}
				}

			} catch (Exception e) {
				//				System.err.println(scopId);
				//				e.printStackTrace();
			} finally {
				i++;
			}
		}

	}

	private void writeToFile() {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(output)));
			pw.print(this);
		} catch (IOException e) {
			throw new RuntimeException("Couldn't output results", e);
		} finally {
			if (pw != null) pw.close();
		}
	}

	public Map<String, String> getFormulas() {
		return formulas;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<String,String> entry : formulas.entrySet()) {
			sb.append(entry.getKey() + "\t" + entry.getValue() + StatUtils.NEWLINE);
		}
		return sb.toString();
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 3) {
			System.err.println("Usage: " + LigandFinder.class.getSimpleName() + " census-file.xml output-file [radius]");
			return;
		}
		int radius = 10;
		if (args.length > 2) {
			radius = Integer.parseInt(args[2]);
		}
		LigandFinder finder = new LigandFinder(radius);
		finder.setOutput(new File(args[1]));
		Results census = Results.fromXML(new File(args[0]));
		finder.find(census);
		System.out.println(finder);
	}

}
