package writers;

import java.io.IOException;

import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;

/**
 * Writes the CeSymm result multiple alignment in FASTA format in a single
 * file. Different entries are split by the characters '//'.
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmFastaWriter extends CeSymmWriter {
	public CeSymmFastaWriter(String filename) throws IOException {
		super(filename);
	}

	@Override
	public synchronized void writeResult(CeSymmResult result) {
		if (result != null ) {
			MultipleAlignment alignment = result.getMultipleAlignment();
			if(alignment != null) {
				writer.write(MultipleAlignmentWriter.toFASTA(alignment));
			}
		}
		writer.println("//");
		writer.flush();
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		// No header for Fasta files
	}
}