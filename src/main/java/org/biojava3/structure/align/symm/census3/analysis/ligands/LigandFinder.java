package org.biojava3.structure.align.symm.census3.analysis.ligands;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomPositionMap;
import org.biojava.bio.structure.AtomPositionMap.GroupMatcher;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava3.structure.align.symm.census3.CensusAlignment;
import org.biojava3.structure.align.symm.census3.CensusResult;
import org.biojava3.structure.align.symm.census3.CensusResultList;
import org.biojava3.structure.align.symm.census3.CensusSignificance;
import org.biojava3.structure.align.symm.census3.CensusSignificanceFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Find ligands near the centroids of symmetric domains.
 * 
 * @author dmyersturnbull
 */
public class LigandFinder {

	private final static Logger logger = LoggerFactory.getLogger(LigandFinder.class);

	private GroupMatcher exclusionMatcher = new GroupMatcher() {
		@Override
		public boolean matches(Group group) {
			ResidueType type = group.getChemComp().getResidueType();
			return group.hasAtom(StructureTools.CA_ATOM_NAME)
					|| AtomPositionMap.AMINO_ACID_NAMES.contains(group.getPDBName())
					|| type == ResidueType.lPeptideLinking || type == ResidueType.glycine
					|| type == ResidueType.lPeptideAminoTerminus || type == ResidueType.lPeptideCarboxyTerminus
					|| type == ResidueType.dPeptideLinking || type == ResidueType.dPeptideAminoTerminus
					|| type == ResidueType.dPeptideCarboxyTerminus;
		}
	};
	private File output;
	private int printFrequency = 100;
	private CensusSignificance significance = CensusSignificanceFactory.forCeSymmOrd();

	public void setSignificance(CensusSignificance significance) {
		this.significance = significance;
	}

	public void setPrintFrequency(int printFrequency) {
		this.printFrequency = printFrequency;
	}

	private LigandList ligandList;

	public LigandFinder() {
	}

	public void setExclusionMatcher(GroupMatcher exclusionMatcher) {
		this.exclusionMatcher = exclusionMatcher;
	}

	public void setOutput(File output) {
		this.output = output;
	}

	private static Atom calcCentroidFromAfpChain(AFPChain afpChain, Atom[] ca) throws StructureException {
		Map<Integer,Integer> map = AlignmentTools.alignmentAsMap(afpChain);
		Atom[] alignedAtoms = new Atom[map.size()];
		int j = 0;
		for (int x : map.keySet()) {
			alignedAtoms[j] = ca[x];
			j++;
		}
		return Calc.getCentroid(alignedAtoms);
	}

	public void find(CensusResultList census) {

		// do this to get a better distribution before we've finished
		Collections.shuffle(census.getEntries());

		/*
		 * include all previous results in the file;
		 * don't redo them
		 */
		if (output != null && output.exists()) {
			try {
				ligandList = LigandList.fromXml(output);
			} catch (IOException e) {
				throw new RuntimeException("Couldn't load " + output);
			}
		} else {
			ligandList = new LigandList();
		}

		AtomCache cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);

		int i = 0;
		for (CensusResult result : census.getEntries()) {

			String scopId = result.getId();
			if (ligandList.contains(scopId)) continue; // don't redo domains

			try {

				/*
				 * Don't even output these.
				 */
				if (!significance.isSignificant(result)) {
					continue;
				}

				/*
				 * Get the structure.
				 */
				Structure structure = cache.getStructure(scopId);
				Atom[] ca = cache.getAtoms(result.getId());

				/*
				 * Find the axis and centroid of the aligned residues
				 */
				Atom centroid;
				RotationAxis axis;
				CensusAlignment mapping = result.getAlignment();
				if (result.getAxis() == null) {
					logger.error("Alignment mapping does not exist");
				}
				try {
					AFPChain afpChain = mapping.buildAfpChain(ca, StructureTools.getAtomCAArray(structure));
					axis = new RotationAxis(afpChain);
					centroid = calcCentroidFromAfpChain(afpChain, ca);
				} catch (Exception e) {
					logger.error("Couldn't use alignment mapping to reconstruct AFPChain", e);
					continue;
				}

				/*
				 *  We want to get all groups near the center that are NOT amino acids or water atoms.
				 *  These are all within the radius of the center.
				 */
				AtomPositionMap atomPositions = new AtomPositionMap(ca, exclusionMatcher);
				Set<ResidueNumber> excluded = atomPositions.getNavMap().keySet();
				Map<Group, Double> distancesToCentroid = StructureTools.getGroupDistancesWithinShell(structure, centroid,
						excluded, Double.POSITIVE_INFINITY, false, false);


				/*
				 * Calculate the distance of every group to the axis.
				 */
				Map<Group,Double> distancesToAxis = new HashMap<Group,Double>();
				for (Group group : distancesToCentroid.keySet()) {
					double minDistance = Double.POSITIVE_INFINITY;
					for (Atom atom : group.getAtoms()) {
						double distance = axis.getProjectedDistance(atom);
						if (distance < minDistance) {
							minDistance = distance;
						}
					}
					distancesToAxis.put(group, minDistance);
				}

				/*
				 *  Add all the ligands to list.
				 */
				if (!distancesToCentroid.isEmpty()) {
					logger.debug("Assigning " + distancesToCentroid.size() + " ligands to " + scopId);
					for (Group group : distancesToCentroid.keySet()) { // constrained by radius
						CensusLigand ligand = new CensusLigand(group.getAtoms(), distancesToCentroid.get(group), distancesToAxis.get(group));
						ligandList.add(scopId, ligand);
					}
				} else {
					logger.debug("Found no center ligands for " + scopId);
					ligandList.put(scopId, new LigandsOfStructure(scopId));
				}

			} catch (Exception e) {
				logger.error(e.getClass().getSimpleName() + " on " + scopId, e);
			} finally {
				i++;
				/*
				 * Print all the results to the file periodically.
				 */
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
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2) {
			System.err.println("Usage: " + LigandFinder.class.getSimpleName() + " census-file.xml output-file");
			return;
		}
		LigandFinder finder = new LigandFinder();
		finder.setOutput(new File(args[1]));
		CensusResultList census = CensusResultList.fromXML(new File(args[0]));
		finder.find(census);
	}

}
