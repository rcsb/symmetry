package writers;

import java.io.IOException;

import org.biojava.nbio.structure.StructureIdentifier;

/**
 * Writer for the Combined Symmetry result.
 * 
 * @author Aleix Lafita
 * 
 */
public class CombinedSymmWriter extends OutputWriter {

	public CombinedSymmWriter(String filename) throws IOException {
		super(filename);
	}

	/**
	 * Writes a line to the file with the Combined Symmetry results of an entry.
	 * 
	 * @param id
	 *            structure identifier
	 * @param chains
	 *            number of chains in the assembly
	 * @param qs
	 *            quaternary symmetry of the assembly
	 * @param is
	 *            true if internal symmetry was found in any of the chains
	 * @param qsis
	 *            symmetry of the combined internal and quaternary
	 * @throws IOException
	 */
	public synchronized void writeResult(StructureIdentifier id, int chains,
			String qs, boolean is, String qsis) throws IOException {

		// Calculate all chains of the BU
		writer.println(id.getIdentifier() + "\t" + chains
				+ "\t" + qs + "\t" + is + "\t" + qsis);
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		writer.println("PDB\t" + "Nr.Chains\t" + "QS\t" + "IS\t" + "QS+IS\t");
		writer.flush();
	}

}