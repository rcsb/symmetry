package org.biojava3.structure.align.symm.census2.analysis;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomPositionMap;
import org.biojava.bio.structure.AtomPositionMap.GroupMatcher;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.model.AFPChain;
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
 * Finds ligands near symmetric interfaces.
 * @author dmyersturnbull
 */
public class InterfaceLigandFinder {

	private static final double DEFAULT_MAX_DISTANCE = 4;

	private static final Logger logger = LogManager.getLogger(InterfaceLigandFinder.class.getName());

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 3) {
			System.err.println("Usage: " + InterfaceLigandFinder.class.getSimpleName()
					+ " census-file.xml output-file [max-distance]");
			return;
		}
		int radius = LigandFinder.DEFAULT_RADIUS;
		if (args.length > 2) {
			radius = Integer.parseInt(args[2]);
		}
		InterfaceLigandFinder finder = new InterfaceLigandFinder();
		finder.setMaxDistance(radius);
		finder.setOutput(new File(args[1]));
		Results census = Results.fromXML(new File(args[0]));
		finder.find(census);
		System.out.println(finder);
	}

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

	private LigandList ligandList;

	private double maxDistance = DEFAULT_MAX_DISTANCE;

	private File output;
	
	private Significance significance;

	private int printFrequency = 100;

	public void setSignificance(Significance significance) {
		this.significance = significance;
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

		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
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
				Structure structure;
				StructureName theName = new StructureName(scopId);
				if (theName.isScopName()) {
					structure = cache.getStructureForDomain(scopId, scop);
				} else {
					structure = cache.getStructure(scopId);
				}
				Atom[] allAtoms = StructureTools.getAllAtomArray(structure);

				// run CE-Symm to get alignment
				FullInfo info = CensusJob.runOn(scopId, Census.AlgorithmGiver.getDefault(), significance, cache, scop);
				AFPChain afpChain = info.getAfpChain();

				/*
				 * add ligands to list
				 */
				StructureLigands ligands = findAlongInterface(scopId, afpChain, allAtoms);
				ligandList.put(scopId, ligands);
				if (!ligands.isEmpty()) {
					logger.debug("Assigning " + ligands.size() + " ligands to " + scopId);
				} else {
					logger.debug("Found no interface ligands for " + scopId);
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

	private StructureLigands findAlongInterface(String scopId, AFPChain afpChain, Atom[] ca) throws StructureException {
		RotationAxis axis = new RotationAxis(afpChain);
		StructureLigands ligands = new StructureLigands(scopId);
		HashSet<Group> groupsFound = new HashSet<Group>();
		for (Atom atom : ca) {
			// only heteroatoms
			if (!exclusionMatcher.matches(atom.getGroup()) && !groupsFound.contains(atom.getGroup())) {
				double distance = axis.getProjectedDistance(atom);
				if (distance <= maxDistance) {
					ligands.add(new Ligand(atom.getGroup().getAtoms(), distance));
					groupsFound.add(atom.getGroup());
				}
			}
		}
		return ligands;
	}

	public Map<String, String> getFormulas() {
		Map<String, String> formulas = new HashMap<String, String>();
		for (StructureLigands ligands : ligandList.values()) {
			formulas.put(ligands.getStructureName(), ligands.toString());
		}
		return formulas;
	}

	public void setExclusionMatcher(GroupMatcher exclusionMatcher) {
		this.exclusionMatcher = exclusionMatcher;
	}

	public void setMaxDistance(double maxDistance) {
		this.maxDistance = maxDistance;
	}

	public void setOutput(File output) {
		this.output = output;
	}

	public void setPrintFrequency(int printFrequency) {
		this.printFrequency = printFrequency;
	}

}
