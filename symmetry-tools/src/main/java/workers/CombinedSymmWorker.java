package workers;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import writers.CombinedSymmWriter;

/**
 * This worker downloads the biological assembly 1 from the PDB, runs the
 * {@link QuatSymmetryDetector} algorithm first and then {@link CeSymm} on each
 * chain of the assembly. Chains are split according to their internal symmetry
 * and the {@link QuatSymmetryDetector} is run again on the complex.
 * <p>
 * The following outcomes are possible: (1) The protein was monomeric and
 * internally symmetric, (2) the complex was heteromeric and asymmetric and
 * gained overall symmetry with internal symmetry, uneven stoichiometry, and (3)
 * the complex was homomeric and symmetric and internal symmetry amplified the
 * point group order.
 * 
 * @author Aleix Lafita
 * @since Apr 2016
 * @version beta
 *
 */
public class CombinedSymmWorker implements Runnable {

	private static final Logger logger = LoggerFactory
			.getLogger(CombinedSymmWorker.class);

	private StructureIdentifier id;
	private CESymmParameters iParams;
	private QuatSymmetryParameters qParams;
	private List<CombinedSymmWriter> writers;

	/**
	 * 
	 * @param id
	 * @param iParams
	 * @param qParams
	 * @param writers
	 */
	public CombinedSymmWorker(StructureIdentifier id, CESymmParameters iParams,
			QuatSymmetryParameters qParams, List<CombinedSymmWriter> writers) {
		this.id = id;
		this.writers = writers;
		this.iParams = iParams;
		this.qParams = qParams;
	}

	@Override
	public void run() {

		try {
			// Obtain the biological assembly 1 of the Structure
			Structure structure;
			try {
				structure = StructureIO.getBiologicalAssembly(id
						.getIdentifier());
			} catch (IOException | StructureException e) {
				logger.error("Could not load Structure " + id.getIdentifier(),
						e);
				return;
			}

			// If there is no biological assembly return
			if (structure == null) {
				logger.warn("Structure " + id.getIdentifier()
						+ " does not have annotated biological assembly.");
				for (CombinedSymmWriter writer : writers)
					writer.writeResult(id, 0, "NA", false, "NA");
				return;
			}

			// Run the Quaternary Symmetry detection (take the optimal)
			QuatSymmetryResults qs = new QuatSymmetryDetector(structure,
					qParams).getGlobalSymmetry().get(0);

			boolean is = false;
			int chains = 0;

			// Run the internal symmetry analysis for every chain and split them
			for (int m = 0; m < structure.nrModels(); m++) {
				List<Chain> newModel = new ArrayList<Chain>(structure.getModel(
						m).size());
				for (int c = 0; c < structure.getModel(m).size(); c++) {
					Atom[] atoms = StructureTools
							.getRepresentativeAtomArray(structure
									.getChain(m, c));
					CeSymmResult result = CeSymm.analyze(atoms, iParams);
					if (result.isSignificant()) {
						is = true;
						Structure divided = SymmetryTools
								.getQuaternaryStructure(result);
						newModel.addAll(divided.getChains());
					} else {
						newModel.add(structure.getChain(m, c));
					}
					chains++;
				}
				structure.setModel(m, newModel);
			}

			// Run the Quaternary Symmetry detection (take the optimal)
			QuatSymmetryResults qsis = new QuatSymmetryDetector(structure,
					qParams).getGlobalSymmetry().get(0);

			// Write into the output files
			for (CombinedSymmWriter writer : writers) {
				try {
					writer.writeResult(id, chains, qs.getSymmetry(), is,
							qsis.getSymmetry());
				} catch (Exception e) {
					logger.error(
							"Could not save results for " + id.getIdentifier(),
							e);
				}
			}

		} catch (Exception e) {
			logger.error("Could not complete job: " + id.getIdentifier(), e);
		} finally {
			logger.info("Finished job: " + id);
		}
	}
}
