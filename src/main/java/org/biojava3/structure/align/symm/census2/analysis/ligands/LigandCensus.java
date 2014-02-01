package org.biojava3.structure.align.symm.census2.analysis.ligands;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomPositionMap;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.AtomPositionMap.GroupMatcher;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Census;
import org.biojava3.structure.align.symm.census2.RandomCensus;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Results;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;

/**
 * A Census that runs only on symmetric ligand-containing domains.
 * @author dmyersturnbull
 */
public class LigandCensus extends Census {

	private static final Logger logger = LogManager.getLogger(LigandCensus.class.getPackage().getName());
	private List<ScopDomain> list;

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
	public static void buildDefault(File censusFile, File inputCensus, boolean shuffle) {
		try {
			int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
			LigandCensus census = new LigandCensus(maxThreads, inputCensus, SignificanceFactory.rotationallySymmetricSmart(), shuffle);
			census.setOutputWriter(censusFile);
			census.setCache(new AtomCache());
			census.run();
			System.out.println(census);
		} catch (Exception e) {
			logger.fatal(e.getMessage(), e);
		}
	}

	public LigandCensus(int maxThreads, File inputCensusFile, Significance sig, boolean shuffle) throws IOException, StructureException {
		super(maxThreads);
		AtomCache cache = new AtomCache();
		Results inputCensus = Results.fromXML(inputCensusFile);
		list = new ArrayList<ScopDomain>();
		for (Result result : inputCensus.getData()) {
			try {
			if (sig.isSignificant(result)) {
				Structure structure = cache.getStructure(result.getScopId());
				Atom[] ca = StructureTools.getAllAtomArray(structure);
				AtomPositionMap atomPositions = new AtomPositionMap(ca, exclusionMatcher);
				Set<ResidueNumber> excluded = atomPositions.getNavMap().keySet();
				Atom centroid = Calc.getCentroid(ca);
				Map<Group, Double> distancesToCentroid = StructureTools.getGroupDistancesWithinShell(structure, centroid,
						excluded, Double.POSITIVE_INFINITY, false, false);
				if (!distancesToCentroid.isEmpty()) {
					System.out.println(result.getScopId());
					list.add(ScopFactory.getSCOP().getDomainByScopID(result.getScopId()));
				}
			}
			} catch (Exception e) {
				logger.error(e);
			}
		}
		if (shuffle) {
			Collections.shuffle(list);
		}
	}

	public static void main(String[] args) {
		if (args.length < 2 || args.length > 3) {
			System.err.println("Usage: " + RandomCensus.class.getSimpleName() + " input-census-file output-census-file [\"shuffle\"]");
			return;
		}
		ScopFactory.setScopDatabase(ScopFactory.getSCOP(ScopFactory.VERSION_1_75A));
		final File inputFile = new File(args[0]);
		final File censusFile = new File(args[1]);
		boolean shuffle = false;
		if (args.length > 2) {
			if (args[2].toLowerCase().equals("shuffle") || args[2].toLowerCase().equals("true")) shuffle = true;
		}
		buildDefault(censusFile, inputFile, shuffle);
	}


}
