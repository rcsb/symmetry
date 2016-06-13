package writers;

import java.io.IOException;

import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;

/**
 * Writes the CeSymm result multiple alignment as aligned residue groups by
 * columns in a single TSV file. This is an easy parsable file of the multiple
 * alignment for other algorithms. Different entries are split by the characters
 * '//'.
 * 
 * @author Aleix Lafita
 * @author Spencer Bliven
 *
 */
public class CeSymmTsvWriter extends CeSymmWriter {

	public CeSymmTsvWriter(String filename) throws IOException {
		super(filename);
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		// no header
	}

	@Override
	public synchronized void writeResult(CeSymmResult result)
			throws IOException {
		if (result != null) {
			MultipleAlignment alignment = result.getMultipleAlignment();
			if (alignment != null)
				writer.write(MultipleAlignmentWriter
						.toAlignedResidues(alignment));
			else {
				// No alignment; just write header
				writer.format("#Struct1:\t%s%n", result.getStructureId());
				writer.format("#Insignificant Alignment%n");
			}
		}
		writer.println("//");
		writer.flush();
	}
}