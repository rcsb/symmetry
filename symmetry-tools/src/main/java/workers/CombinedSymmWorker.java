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
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.symmetry.core.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorPointGroup;
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
	private boolean show3d;

	/**
	 * 
	 * @param id
	 * @param iParams
	 * @param qParams
	 * @param writers
	 * @param show3d
	 */
	public CombinedSymmWorker(StructureIdentifier id, CESymmParameters iParams,
			QuatSymmetryParameters qParams, List<CombinedSymmWriter> writers,
			boolean show3d) {
		this.id = id;
		this.writers = writers;
		this.iParams = iParams;
		this.qParams = qParams;
		this.show3d = show3d;
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
				for (CombinedSymmWriter writer : writers)
					writer.writeBlank(id);
				return;
			}

			// If there is no biological assembly return
			if (structure == null) {
				logger.error("Structure " + id.getIdentifier()
						+ " does not have annotated biological assembly.");
				for (CombinedSymmWriter writer : writers)
					writer.writeBlank(id);
				return;
			}
			structure.setStructureIdentifier(id);

			// Run the Quaternary Symmetry detection (take the optimal)
			QuatSymmetryResults qs = new QuatSymmetryDetector(structure,
					qParams).getGlobalSymmetry().get(0);

			boolean is = false;
			int chainNr = 0;
			char chainId = 'A';

			// Run the internal symmetry analysis for every chain and split them
			for (int m = 0; m < structure.nrModels(); m++) {
				List<Chain> newModel = new ArrayList<Chain>(structure.getModel(
						m).size());
				for (int c = 0; c < structure.getModel(m).size(); c++) {
					Atom[] atoms = StructureTools
							.getRepresentativeAtomArray(structure
									.getChainByIndex(m, c));
					CeSymmResult result = CeSymm.analyze(atoms, iParams);
					if (result.isSignificant()) {
						is = true;
						List<Structure> repeats = SymmetryTools
								.divideStructure(result);
						for (Structure s : repeats) {
							newModel.add(s.getChainByIndex(0));
						}
					} else {
						structure.getChainByIndex(m, c).setChainID((chainId++) + "");
						newModel.add(structure.getChainByIndex(m, c));
					}
					chainNr++;
				}
				structure.setModel(m, newModel);
			}

			// Run the Quaternary Symmetry detection (take the pseudosymmetric)
			List<QuatSymmetryResults> qsiss = new QuatSymmetryDetector(
					structure, qParams).getGlobalSymmetry();
			QuatSymmetryResults qsis = qsiss.get(0);
			for (QuatSymmetryResults q : qsiss) {
				if (q.getSubunits().isPseudoSymmetric())
					qsis = q;
			}

			if (show3d) {
				AxisAligner aligner = AxisAligner.getInstance(qsis);
				JmolSymmetryScriptGenerator scriptGenerator = JmolSymmetryScriptGeneratorPointGroup
						.getInstance(aligner, "g");
				String script = "set defaultStructureDSSP true; set measurementUnits ANGSTROMS;  select all;  spacefill off; wireframe off; "
						+ "backbone off; cartoon on; color cartoon structure; color structure;  select ligand;wireframe 0.16;spacefill 0.5; "
						+ "color cpk ; select all; model 0;set antialiasDisplay true; autobond=false;save STATE state_1;";
				script += scriptGenerator.getOrientationWithZoom(0);
				script += scriptGenerator.drawPolyhedron();
				script += scriptGenerator.drawAxes();
				script += scriptGenerator.colorBySymmetry();
				script += "draw axes* on; draw poly* on;";

				StructureAlignmentJmol jmol = new StructureAlignmentJmol();
				jmol.setStructure(structure);
				jmol.evalString(script);
			}

			// Write into the output files
			for (CombinedSymmWriter writer : writers) {
				try {
					writer.writeResult(id, chainNr, qs, is, qsis);
				} catch (Exception e) {
					logger.error(
							"Could not save results for " + id.getIdentifier(),
							e);
					writer.writeBlank(id);
				}
			}

		} catch (Exception e) {
			logger.error("Could not complete job: " + id.getIdentifier(), e);
			for (CombinedSymmWriter writer : writers) {
				writer.writeBlank(id);
			}
		} finally {
			logger.info("Finished job: " + id);
		}
	}
}
