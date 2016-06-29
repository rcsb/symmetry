package workers;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.gui.BiojavaJmol;
import org.biojava.nbio.structure.symmetry.axis.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorPointGroup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import writers.QuatSymmWriter;

/**
 * This Runnable implementation runs the Quaternary Symmetry Detector on the
 * input structure and with the input parameters and writes the results to the
 * output writers provided.
 * <p>
 * If the 3D visualization is turned on, it creates a new thread with the Jmol
 * frame.
 * 
 * @author Aleix Lafita
 *
 */
public class QuatSymmWorker implements Runnable {

	private static final Logger logger = LoggerFactory
			.getLogger(QuatSymmWorker.class);

	private StructureIdentifier id;
	private SubunitClustererParameters cparams;
	private QuatSymmetryParameters sparams;
	private AtomCache cache;
	private List<QuatSymmWriter> writers;
	private boolean show3d;

	public QuatSymmWorker(StructureIdentifier id,
			QuatSymmetryParameters sparams, SubunitClustererParameters cparams,
			AtomCache cache, List<QuatSymmWriter> writers, boolean show3d) {
		this.id = id;
		this.cache = cache;
		this.writers = writers;
		this.sparams = sparams;
		this.cparams = cparams;
		this.show3d = show3d;
	}

	@Override
	public void run() {

		try {
			// Obtain the structure representation
			Structure structure = null;
			try {
				structure = cache.getStructure(id);
			} catch (IOException | StructureException e) {
				logger.error("Could not load Structure " + id.getIdentifier(),
						e);
				return;
			}

			// Calculate the global symmetry
			QuatSymmetryResults result = QuatSymmetryDetector
					.calcGlobalSymmetry(structure, sparams, cparams);

			if (result.getSymmetry().equals("C1")) {
				// Calculate local symmetry
				List<QuatSymmetryResults> localResults = QuatSymmetryDetector
						.calcLocalSymmetries(structure, sparams, cparams);
				QuatSymmetryResults local = null;
				for (QuatSymmetryResults r : localResults) {
					if (local == null)
						local = r;
					else if (local.getSubunits().getSubunitCount() < r
							.getSubunits().getSubunitCount())
						local = r;
				}
				if (local != null)
					result = local;
			}
			
			// Write into the output files
			for (QuatSymmWriter writer : writers) {
				try {
					synchronized (writer) {
						writer.writeResult(id.toString(), result);
					}
				} catch (Exception e) {
					logger.error(
							"Could not save results for " + id.getIdentifier(),
							e);
				}
			}

			if (show3d && result != null) {
				AxisAligner aligner = AxisAligner.getInstance(result);
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

				BiojavaJmol jmol = new BiojavaJmol();
				jmol.setStructure(structure);
				jmol.evalString(script);
			}

		} catch (Exception e) {
			logger.error("Could not complete job: " + id.getIdentifier(), e);
		} finally {
			logger.info("Finished job: " + id);
		}
	}
}
