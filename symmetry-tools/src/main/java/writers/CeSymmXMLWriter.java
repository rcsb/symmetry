package writers;

import java.io.IOException;

import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;

/**
 * Writes the XML serialization of the multiple structure alignment of all
 * entries in a single file.
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmXMLWriter extends CeSymmWriter {

	public CeSymmXMLWriter(String filename) throws IOException {
		super(filename);
	}

	@Override
	public synchronized void writeResult(CeSymmResult result)
			throws IOException {
		if (result != null && result.getMultipleAlignment() != null) {
			writer.append(MultipleAlignmentWriter.toXML(result
					.getMultipleAlignment().getEnsemble()));
			writer.flush();
		}
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		// No header for XML file
	}

}