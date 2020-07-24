package writers;

import java.io.IOException;

import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;

/**
 * Parent class for all QuatSymm result output formats.
 * 
 * @author Aleix Lafita
 * 
 */
public abstract class QuatSymmWriter extends OutputWriter {

	public QuatSymmWriter(String filename) throws IOException {
		super(filename);
	}

	/**
	 * Writes a line to the file with the QuatSymm results of an entry.
	 * Implementations of this method need to be synchronized to avoid writting
	 * at the same time.
	 * 
	 * @param identifier
	 * @param result
	 * @throws IOException
	 */
	abstract public void writeResult(String identifier,
			QuatSymmetryResults result) throws Exception;

}