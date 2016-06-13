package writers;

import java.io.IOException;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;

/**
 * Parent class for all CeSymm result output formats.
 * 
 * @author Aleix Lafita
 * 
 */
public abstract class CeSymmWriter extends OutputWriter {

	public CeSymmWriter(String filename) throws IOException {
		super(filename);
	}

	/**
	 * Writes a line to the file with the CeSymm results of an entry.
	 * Implementations of this method need to be synchronized to avoid writting
	 * at the same time.
	 * 
	 * @param result
	 * @throws IOException
	 */
	abstract public void writeResult(CeSymmResult result) throws IOException;

}