package writers;

import java.io.IOException;

import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;

/**
 * Writes the CeSymm result multiple alignment in FatCat format in a single
 * file. Different entries are split by the characters '//'.
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmFatcatWriter extends CeSymmWriter {

	public CeSymmFatcatWriter(String filename) throws IOException {
		super(filename);
	}

	@Override
	public synchronized void writeResult(CeSymmResult result) {
		if (result != null) {
			MultipleAlignment alignment = result.getMultipleAlignment();
			if (alignment != null) {
				writer.write(MultipleAlignmentWriter.toFatCat(alignment));
			} else {
				writer.format("Structures:[%s]%n", result.getStructureId());
				writer.format("Insignificant Alignment%n");
			}
		}
		writer.println("//");
		writer.flush();
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		// No header for FatCat file
	}
}