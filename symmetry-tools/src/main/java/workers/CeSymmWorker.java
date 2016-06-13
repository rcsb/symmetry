package workers;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.symmetry.gui.SymmetryDisplay;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import writers.CeSymmWriter;

/**
 * This Runnable implementation runs CeSymm on the input structure and with the
 * input parameters and writes the results to the output writers provided.
 * <p>
 * If the 3D visualization is turned on, it creates a new thread with the Jmol
 * frame.
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmWorker implements Runnable {
	
	private static final Logger logger = LoggerFactory
			.getLogger(CeSymmWorker.class);

	private StructureIdentifier id;
	private CESymmParameters params;
	private AtomCache cache;
	private List<CeSymmWriter> writers;
	private boolean show3d;

	public CeSymmWorker(StructureIdentifier id, CESymmParameters params,
			AtomCache cache, List<CeSymmWriter> writers, boolean show3d) {
		this.id = id;
		this.cache = cache;
		this.writers = writers;
		this.show3d = show3d;
		this.params = params;
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
			
			Atom[] atoms = SymmetryTools.getRepresentativeAtoms(structure);
			
			// Run the symmetry analysis
			CeSymmResult result = CeSymm.analyze(atoms, params);

			// Write into the output files
			for (CeSymmWriter writer : writers) {
				try {
					synchronized (writer) {
						writer.writeResult(result);
					}
				} catch (Exception e) {
					logger.error(
							"Could not save results for " + id.getIdentifier(),
							e);
				}
			}

			// Display alignment in 3D Jmol
			if (show3d) {
				SymmetryDisplay.display(result);
			}
		} catch (Exception e) {
			logger.error("Could not complete job: " + id.getIdentifier(), e);
		} finally {
			logger.info("Finished job: " + id);
		}
	}
}