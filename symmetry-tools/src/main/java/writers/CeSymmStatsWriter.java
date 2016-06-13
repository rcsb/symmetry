package writers;

import java.io.IOException;

import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;

/**
 * Writes a stats summary of the CeSymm result in TSV format.
 * 
 * @author Aleix Lafita
 * @author Spencer Bliven
 *
 */
public class CeSymmStatsWriter extends CeSymmWriter {

	public CeSymmStatsWriter(String filename) throws IOException {
		super(filename);
	}

	@Override
	public synchronized void writeHeader() {
		writer.println("Name\t" + "NumRepeats\t" + "SymmGroup\t" + "Refined\t"
				+ "SymmLevels\t" + "SymmType\t" + "RotationAngle\t"
				+ "ScrewTranslation\t" + "UnrefinedTMscore\t"
				+ "UnrefinedRMSD\t" + "SymmTMscore\t" + "SymmRMSD\t"
				+ "RepeatLength\t" + "CoreLength\t" + "Length\t" + "Coverage");
		writer.flush();
	}

	@Override
	public synchronized void writeResult(CeSymmResult result)
			throws IOException {
		String id = null;
		if (result == null) {
			writeEmptyRow(id);
			writer.flush();
			return;
		}
		try {
			id = result.getStructureId().getIdentifier();
			int repeatLen = 0;
			int coreLen = result.getSelfAlignment().getOptLength();
			double coverage = result.getSelfAlignment().getCoverage1() / 100;
			String group = result.getSymmGroup();
			int order = result.getNumRepeats();
			double symmrmsd = 0.0;
			double symmscore = 0.0;

			RotationAxis rot = new RotationAxis(result.getSelfAlignment()
					.getBlockRotationMatrix()[0], result.getSelfAlignment()
					.getBlockShiftVector()[0]);
			double rotation_angle = Math.toDegrees(rot.getAngle());
			double screw_translation = new Vector3d(rot.getScrewTranslation()
					.getCoords()).length();

			int structureLen = result.getAtoms().length;

			// If there is refinement alignment
			if (result.isRefined()) {
				MultipleAlignment msa = result.getMultipleAlignment();
				symmrmsd = msa.getScore(MultipleAlignmentScorer.RMSD);
				symmscore = msa.getScore(MultipleAlignmentScorer.AVGTM_SCORE);

				repeatLen = msa.length();
				coreLen = msa.getCoreLength() * msa.size();

				// Calculate coverage
				coverage = 0;
				for (int s = 0; s < msa.size(); s++)
					coverage += msa.getCoverages().get(s);

				rot = new RotationAxis(result.getAxes().getElementaryAxes()
						.get(0));
				rotation_angle = Math.toDegrees(rot.getAngle());
				screw_translation = new Vector3d(rot.getScrewTranslation()
						.getCoords()).length();
			}

			writer.format("%s\t%d\t%s\t%b\t%d\t%s\t%.2f\t%.2f\t%.2f\t"
					+ "%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f%n", id, order, group,
					result.isRefined(), result.getSymmLevels(), result
							.getAxes().getElementaryAxis(0).getSymmType(),
					rotation_angle, screw_translation, result
							.getSelfAlignment().getTMScore(), result
							.getSelfAlignment().getTotalRmsdOpt(), symmscore,
					symmrmsd, repeatLen, coreLen, structureLen, coverage);
		} catch (Exception e) {
			// If any exception occurs when writing the results store empty row
			logger.warn("Could not write result for entry: " + id
					+ ". Writting empty row.");
			writeEmptyRow(id);
		}

		writer.flush();
	}

	private synchronized void writeEmptyRow(String id) {
		writer.format("%s\t%d\t%s\t%b\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t"
				+ "%.2f\t%d\t%d\t%d\t%.2f%n", id, 1, "C1", false, 0,
				SymmetryType.DEFAULT, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0,
				0.0);
	}
}