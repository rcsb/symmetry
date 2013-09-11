package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

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

/**
 * Find ligands near the centroids of symmetric domains.
 * 
 * @author dmyerstu
 */
public class LigandFinder {

	public static int DEFAULT_RADIUS = 5;

	private static final Logger logger = LogManager.getLogger(LigandFinder.class.getName());

	private GroupMatcher exclusionMatcher = new GroupMatcher() {
		@Override
		public boolean matches(Group group) {
			ResidueType type = group.getChemComp().getResidueType();
			return group.hasAtom(StructureTools.caAtomName)
					|| AtomPositionMap.AMINO_ACID_NAMES.contains(group.getPDBName())
					|| type == ResidueType.lPeptideLinking || type == ResidueType.glycine
					|| type == ResidueType.lPeptideAminoTerminus || type == ResidueType.lPeptideCarboxyTerminus
					|| type == ResidueType.dPeptideLinking || type == ResidueType.dPeptideAminoTerminus
					|| type == ResidueType.dPeptideCarboxyTerminus;
		}
	};
	private int radius = DEFAULT_RADIUS;
	private File output;
	private int printFrequency = 20;

	public void setPrintFrequency(int printFrequency) {
		this.printFrequency = printFrequency;
	}

	private LigandList ligandList;

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

	public void find(Results census) {

		if (output != null && output.exists()) {
			try {
				ligandList = LigandList.fromXml(output);
			} catch (IOException e) {
				throw new RuntimeException("Couldn't load " + output);
			}
		} else {
			ligandList = new LigandList();
		}

		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
		Significance sig = SignificanceFactory.rotationallySymmetricSmart();
		AtomCache cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);

		int i = 0;
		for (Result result : census.getData()) {

			String scopId = result.getScopId();
			if (ligandList.contains(scopId)) continue; // don't redo domains
			try {

				ScopDomain domain = scop.getDomainByScopID(scopId);
				if (domain == null) {
					logger.error(result.getScopId() + " is null");
					continue;
				}
				if (!sig.isSignificant(result)) {
					continue;
				}

				// we want to get all groups in the center that are NOT amino
				// acids or water atoms
				Structure structure = cache.getStructureForDomain(scopId, scop);
				Atom[] ca = StructureTools.getAtomCAArray(structure);
				Atom centroid = Calc.getCentroid(ca);
				AtomPositionMap atomPositions = new AtomPositionMap(ca, exclusionMatcher);
				Set<ResidueNumber> excluded = atomPositions.getNavMap().keySet();
				// Set<Group> ligands =
				// StructureTools.getGroupsWithinShell(structure, centroid,
				// residues, radius, false);
				Map<Group, Double> ligandDistances = StructureTools.getGroupDistancesWithinShell(structure, centroid,
						excluded, radius, false, false);

				if (!ligandDistances.isEmpty()) {
					for (Group group : ligandDistances.keySet()) {
						ligandList.put(scopId, new Ligand(group.getAtoms(), ligandDistances.get(group)));
						if (i % printFrequency == 0 && i > 0 && output != null) {
							logger.info("Printing to " + output.getPath());
							try {
								ligandList.writeXmlToFile(output);
							} catch (IOException e1) {
								logger.error("Couldn't write to " + output, e1);
							}
						}
					}
				}

			} catch (Exception e) {
				e.printStackTrace();
				logger.error(e.getClass().getSimpleName() + " on " + scopId);
			} finally {
				i++;
			}
		}

	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 3) {
			System.err
					.println("Usage: " + LigandFinder.class.getSimpleName() + " census-file.xml output-file [radius]");
			return;
		}
		int radius = LigandFinder.DEFAULT_RADIUS;
		if (args.length > 2) {
			radius = Integer.parseInt(args[2]);
		}
		LigandFinder finder = new LigandFinder(radius);
		finder.setOutput(new File(args[1]));
		Results census = Results.fromXML(new File(args[0]));
		finder.find(census);
		System.out.println(finder);
	}

	public Map<String, String> getFormulas() {
		Map<String, String> formulas = new HashMap<String, String>();
		for (StructureLigands ligands : ligandList.values()) {
			formulas.put(ligands.getStructureName(), ligands.toString());
		}
		return formulas;
	}

}
