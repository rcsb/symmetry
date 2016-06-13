package writers;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.internal.SymmetryAxes;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.internal.SymmetryAxes.Axis;

/**
 * Prints the SymmetryAxes information.
 * 
 * @author Spencer Bliven
 *
 */
public class CeSymmAxesWriter extends CeSymmWriter {

	public CeSymmAxesWriter(String filename) throws IOException {
		super(filename);
	}

	@Override
	public void writeHeader() {

		writer.println("Structure\t" + "SymmLevel\t" + "SymmType\t" + "SymmGroup\t"
				+ "RotationAngle\t" + "ScrewTranslation\t" + "Point1\t" + "Point2\t"
				+ "AlignedRepeats");
		writer.flush();
	}
	private void writeEmptyRow(String id) {
		writer.format("%s\t%d\t%s\t%s"
				+ "%.2f\t%.2f\t%.3f,%.3f,%.3f\t%.3f,%.3f,%.3f\t%s%n",
				id, -1, SymmetryType.DEFAULT, "C1",
				0., 0., 0.,0.,0., 0.,0.,0.,
				"");
	}
	@Override
	public void writeResult(CeSymmResult result) throws IOException {
		String id = null;
		if( result == null) {
			writeEmptyRow(id);
			writer.flush();
			return;
		}
		try {
			id = result.getStructureId().getIdentifier();
			
			SymmetryAxes axes = result.getAxes();
			Atom[] atoms = result.getAtoms();
			//TODO is this correct for hierarchical cases?
			String symmGroup = result.getSymmGroup();
			
			for(Axis axis : axes.getSymmetryAxes() ) {
				RotationAxis rot = axis.getRotationAxis();
//				Set<Integer> repIndex = new TreeSet<Integer>(axes
//						.getRepeatRelation(a).get(0));
				List<List<Integer>> repeatsCyclicForm = axes.getRepeatsCyclicForm(axis);
				String cyclicForm = SymmetryAxes.getRepeatsCyclicForm(repeatsCyclicForm, result.getRepeatsID());
				
				// Get atoms aligned by this axis
				List<Atom> axisAtoms = new ArrayList<Atom>();
				for(List<Integer> cycle : repeatsCyclicForm) {
					for(Integer repeat : cycle) {
						axisAtoms.addAll(Arrays.asList(atoms[repeat]));
					}
				}
				Pair<Atom> bounds = rot.getAxisEnds(axisAtoms.toArray(new Atom[axisAtoms.size()]));
						
				Atom start = bounds.getFirst();
				Atom end = bounds.getSecond();

				writer.format("%s\t%d\t%s\t%s\t"
						+ "%.2f\t%.2f\t%.3f,%.3f,%.3f\t%.3f,%.3f,%.3f\t%s%n",
						id, axis.getLevel()+1, axis.getSymmType(), symmGroup,
						Math.toDegrees(rot.getAngle()), rot.getTranslation(),
						start.getX(),start.getY(),start.getZ(),
						end.getX(),end.getY(),end.getZ(),
						cyclicForm);
			}
		} catch (Exception e) {
			// If any exception occurs when writing the results store empty
			// better
			logger.warn("Could not write result... storing empty row.", e);
			writeEmptyRow(id);
		}

		writer.flush();
	}
}