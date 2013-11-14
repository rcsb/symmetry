package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
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
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.align.symm.census2.CensusJob;
import org.biojava3.structure.align.symm.census2.CensusJob.FullInfo;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;

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
	private int printFrequency = 100;
	private boolean useOnlyAligned = true;
	private Significance significance;
	

	public void setSignificance(Significance significance) {
		this.significance = significance;
	}

	public void setUseOnlyAligned(boolean useOnlyAligned) {
		this.useOnlyAligned = useOnlyAligned;
	}

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

		Collections.shuffle(census.getData());
		
		if (output != null && output.exists()) {
			try {
				ligandList = LigandList.fromXml(output);
			} catch (IOException e) {
				throw new RuntimeException("Couldn't load " + output);
			}
		} else {
			ligandList = new LigandList();
		}

		ScopDatabase scop = ScopFactory.getSCOP();
		AtomCache cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);

		int i = 0;
		for (Result result : census.getData()) {

			String scopId = result.getScopId();
			if (ligandList.contains(scopId)) continue; // don't redo domains
			try {

				if (!significance.isSignificant(result)) {
					continue;
				}

				// we want to get all groups in the center that are NOT amino
				// acids or water atoms
//				Structure structure = cache.getStructureForDomain(scopId, scop);
				Structure structure;
				StructureName theName = new StructureName(scopId);
				if (theName.isScopName()) {
					if (scop == null) scop = ScopFactory.getSCOP();
					structure = cache.getStructureForDomain(scopId, scop);
				} else {
					structure = cache.getStructure(scopId);
				}
				Atom[] ca = StructureTools.getAtomCAArray(structure);

				Atom centroid;
				if (useOnlyAligned) {
					// run CE-Symm to get alignment
					FullInfo info = CensusJob.runOn(scopId, Census.AlgorithmGiver.getDefault(), significance, cache, scop);
					RotationAxis axis = new RotationAxis(info.getAfpChain());
					centroid = axis.getRotationPos();
				} else {
					centroid = Calc.getCentroid(ca);
				}

				AtomPositionMap atomPositions = new AtomPositionMap(ca, exclusionMatcher);
				Set<ResidueNumber> excluded = atomPositions.getNavMap().keySet();
				// Set<Group> ligands =
				// StructureTools.getGroupsWithinShell(structure, centroid,
				// residues, radius, false);
				Map<Group, Double> ligandDistances = StructureTools.getGroupDistancesWithinShell(structure, centroid,
						excluded, radius, false, false);

				/*
				 *  add ligands to list
				 */
				if (!ligandDistances.isEmpty()) {
					logger.debug("Assigning " + ligandDistances.size() + " ligands to " + scopId);
					for (Group group : ligandDistances.keySet()) {
						ligandList.add(scopId, new Ligand(group.getAtoms(), ligandDistances.get(group)));
					}
				} else {
					logger.debug("Found no center ligands for " + scopId);
					ligandList.put(scopId, new StructureLigands(scopId));
				}

			} catch (Exception e) {
				e.printStackTrace();
				logger.error(e.getClass().getSimpleName() + " on " + scopId);
			} finally {
				i++;
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

	}

	/**
	 * Removes unaligned groups from {@code ligandDistances}.
	 * @param ligandDistances The map to remove from
	 * @param alignment A map from group position numbers to group position numbers
	 * @param ca All atoms of the structure
	 * @throws StructureException
	 * @throws IOException
	 * @deprecated Irrelevant
	 */
	@Deprecated
	private void removeUnaligned(Map<Group, Double> ligandDistances, Map<Integer,Integer> alignment, Atom[] ca) throws StructureException, IOException {

		// get a list of all of the groups in the domain
		// the index of the list is critical
		List<Group> groupsInStructure = new ArrayList<Group>();
		HashSet<WeakReference<Group>> groupsInStructureSet = new HashSet<WeakReference<Group>>(); // so we don't have to iterate thru in liear time
		for (Atom atom : ca) {
			Group group = atom.getGroup();
			if (!groupsInStructureSet.contains(group)) {
				groupsInStructure.add(group);
				groupsInStructureSet.add(new WeakReference<Group>(group));
			}
		}

		// now find aligned residues
		HashSet<Group> alignedGroups = new HashSet<Group>();
		for (Integer key : alignment.keySet()) {
			alignedGroups.add(groupsInStructure.get(key));
		}

		// now remove any group that's not on the list
		List<Group> toRemove = new ArrayList<Group>();
		for (Group group : ligandDistances.keySet()) {
			if (!alignedGroups.contains(group)) {
				toRemove.add(group);
			}
		}
		for (Group group : toRemove) {
			ligandDistances.remove(group);
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
