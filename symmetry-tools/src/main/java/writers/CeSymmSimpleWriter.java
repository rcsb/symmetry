package writers;

import java.io.IOException;

import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;

/**
 * This is a simple writer meant for the standard out that prints the most
 * relevant CeSymm prediction and a reason for it.
 * 
 * @author Aleix Lafita
 * @author Spencer Bliven
 *
 */
public class CeSymmSimpleWriter extends CeSymmWriter {

	public CeSymmSimpleWriter(String filename) throws IOException {
		super(filename);
	}

	@Override
	public synchronized void writeHeader() {
		writer.println("Structure\tNumRepeats\tSymmGroup\tReason");
		writer.flush();
	}

	private synchronized void writeEmptyRow(String id) {
		writer.format("%s\t%d\t%s\t%s%n", id, 1, "C1", "Error");
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
			writer.append(id);
			writer.append("\t");
			writer.append(Integer.toString(result.getSymmOrder()));
			writer.append("\t");
			writer.append(result.getSymmGroup());
			writer.append("\t");
			writer.append(result.getReason());
			writer.println();
			writer.flush();
		} catch (Exception e) {
			logger.warn("Could not write result for entry: " + id
					+ ". Writting empty row.");
			writeEmptyRow(id);
		}
		writer.flush();
	}
}