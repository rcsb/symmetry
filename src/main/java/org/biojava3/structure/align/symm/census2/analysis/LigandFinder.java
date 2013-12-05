package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomPositionMap;
import org.biojava.bio.structure.AtomPositionMap.GroupMatcher;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.AlignmentMapping;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.align.symm.census2.CensusJob;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;

/**
 * Find ligands near the centroids of symmetric domains.
 * 
 * @author dmyerstu
 */
public class LigandFinder {

	public static double DEFAULT_RADIUS = Double.POSITIVE_INFINITY;

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
	private double radius = DEFAULT_RADIUS;
	private File output;
	private int printFrequency = 100;
	private boolean rebuildMissingAlignments = false;
	private Significance significance;
	

	public void setSignificance(Significance significance) {
		this.significance = significance;
	}

	public void setRebuildMissingAlignments(boolean useOnlyAligned) {
		this.rebuildMissingAlignments = useOnlyAligned;
	}

	public void setPrintFrequency(int printFrequency) {
		this.printFrequency = printFrequency;
	}

	private LigandList ligandList;

	public LigandFinder(double radius) {
		this.radius = radius;
	}

	public LigandFinder(double radius, GroupMatcher exclusionMatcher) {
		this.radius = radius;
		this.exclusionMatcher = exclusionMatcher;
	}

	public void setExclusionMatcher(GroupMatcher exclusionMatcher) {
		this.exclusionMatcher = exclusionMatcher;
	}

	public void setRadius(double radius) {
		this.radius = radius;
	}

	public void setOutput(File output) {
		this.output = output;
	}

	public void find(Results census) {

		// do this to get a better distribution while still running
		Collections.shuffle(census.getData());
		
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

				/*
				 * Get the structure.
				 */
				Structure structure;
				StructureName theName = new StructureName(scopId);
				if (theName.isScopName()) {
					if (scop == null) scop = ScopFactory.getSCOP();
					structure = cache.getStructureForDomain(scopId, scop);
				} else {
					structure = cache.getStructure(scopId);
				}
				Atom[] ca = StructureTools.getAtomCAArray(structure);

				/*
				 * Find the axis and centroid of the aligned residues
				 */
				Atom centroid;
				RotationAxis axis;
				AlignmentMapping mapping = result.getAlignmentMapping();
				if (mapping != null) {
					AFPChain afpChain = mapping.buildAfpChain(StructureTools.getAtomCAArray(structure), StructureTools.getAtomCAArray(structure));
					axis = new RotationAxis(afpChain);
				} else if (rebuildMissingAlignments) {
					// run CE-Symm to get alignment
					CensusJob job = CensusJob.setUpJob(scopId, 0, Census.AlgorithmGiver.getDefault(), significance, cache, scop);
					job.setStoreAfpChain(true);
					job.call();
					axis = new RotationAxis(job.getAfpChain());
				} else {
					logger.warn("Skipping " + scopId + " because the axis could not be found");
					continue;
				}
				centroid = axis.getRotationPos();

				/*
				 *  We want to get all groups near the center that are NOT amino acids or water atoms.
				 *  These are all within the radius of the center.
				 */
				AtomPositionMap atomPositions = new AtomPositionMap(ca, exclusionMatcher);
				Set<ResidueNumber> excluded = atomPositions.getNavMap().keySet();
				Map<Group, Double> distancesToCentroid = StructureTools.getGroupDistancesWithinShell(structure, centroid,
						excluded, radius, false, false);


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
					for (Group group : distancesToCentroid.keySet()) {
						Ligand ligand = new Ligand(group.getAtoms(), distancesToCentroid.get(group), distancesToAxis.get(group));
						ligandList.add(scopId, ligand);
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
		if (args.length < 2 || args.length > 3) {
			System.err
			.println("Usage: " + LigandFinder.class.getSimpleName() + " census-file.xml output-file [radius]");
			return;
		}
		double radius = LigandFinder.DEFAULT_RADIUS;
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
