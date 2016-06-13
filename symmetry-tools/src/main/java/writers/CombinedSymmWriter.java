package writers;

import java.io.IOException;

import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;

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
	 *            quaternary symmetry results of the assembly
	 * @param is
	 *            true if internal symmetry was found in any of the chains
	 * @param qsis
	 *            symmetry results of the combined internal and quaternary
	 * @throws IOException
	 */
	public synchronized void writeResult(StructureIdentifier id, int chains,
			QuatSymmetryResults qs, boolean is, QuatSymmetryResults qsis)
			throws IOException {

		// Calculate all chains of the BU
		writer.println(id.getIdentifier() + "\t" + chains + "\t"
				+ qs.getSymmetry() + "\t" + qs.getSubunits().getStoichiometry()
				+ "\t" + is + "\t" + qsis.getSymmetry() + "\t"
				+ qsis.getSubunits().getStoichiometry());
		writer.flush();
	}
	
	public synchronized void writeBlank(StructureIdentifier id) {
		writer.println(id.getIdentifier() + "\t0\tNA\tNA\tfalse\tNA\tNA");
		writer.flush();
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		writer.println("PDB\t" + "Nr.Ch\t" + "QS-G\t" + "QS-S\t"
				+ "IS?\t" + "IS-G\t" + "IS-S");
		writer.flush();
	}

}